#include <RcppArmadillo.h>
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel, RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;


// Helper Functions
inline double manual_dot(const vec& a, const vec& b) {
  return a(0) * b(0) + a(1) * b(1) + a(2) * b(2);
}

inline vec manual_cross(const vec& a, const vec& b) {
  vec res(3);
  res(0) = a(1) * b(2) - a(2) * b(1);
  res(1) = a(2) * b(0) - a(0) * b(2);
  res(2) = a(0) * b(1) - a(1) * b(0);
  return res;
}

// AABB
struct AABB {
  vec3 min_pt;
  vec3 max_pt;
  int face_idx;
};

// Distance from a point to an AABB (squared)
double dist_sq_point_aabb(const vec3& p, const AABB& box) {
  double sq_dist = 0.0;
  for(int i = 0; i < 3; ++i) {
    if(p[i] < box.min_pt[i]) sq_dist += pow(box.min_pt[i] - p[i], 2);
    else if(p[i] > box.max_pt[i]) sq_dist += pow(p[i] - box.max_pt[i], 2);
  }
  return sq_dist;
}

// Möller-Trumbore algorithm
double rayIntersectsTriangle(const vec& rayOrigin, const vec& rayVector, const vec& v0, const vec& v1, const vec& v2, double epsilon) {
  vec3 edge1 = v1 - v0;
  vec3 edge2 = v2 - v0;

  // backface culling = disabled
  // vec normal = manual_cross(edge1, edge2);
  // if (manual_dot(normal, rayVector) > 0) return -1.0;

  // ray is parallel
  vec3 ray_cross_e2 = manual_cross(rayVector, edge2);
  double det = manual_dot(edge1, ray_cross_e2);
  if (abs(det) < epsilon) return -1.0;

  // ray passes outside edge2's bounds
  double inv_det = 1.0 / det;
  vec3 s = rayOrigin - v0;
  double u = inv_det * manual_dot(s, ray_cross_e2);
  if (u < -epsilon || u > 1.0 + epsilon) return -1.0;

  // ray passes outside edge1's bounds
  vec3 s_cross_e1 = manual_cross(s, edge1);
  double v = inv_det * manual_dot(rayVector, s_cross_e1);
  if (v < -epsilon || u + v > 1.0 + epsilon) return -1.0;

  // if ray line intersects with the triangle
  double t = inv_det * manual_dot(edge2, s_cross_e1);
  return (t > epsilon) ? t : -1.0;
}

// RcppParallel worker for intersection
struct IntersectionWorker : public Worker {
  const mat& centroids;
  const mat& normals;
  const mat& verticesB;
  const imat& facesB;
  const mat& bMins;
  const mat& bMaxs;
  double epsilon;

  // raw pointer to ensure memory is hit
  int* output;

  IntersectionWorker(const mat& c, const mat& n, const mat& vB, const imat& fB,
                     const mat& bMi, const mat& bMa, const double& eps, int* res)
    : centroids(c), normals(n), verticesB(vB), facesB(fB),
      bMins(bMi), bMaxs(bMa), epsilon(eps), output(res) {}

  void operator()(std::size_t begin, std::size_t end) {
    int nB = facesB.n_cols;

    for (std::size_t i = begin; i < end; ++i) {
      vec P = centroids.col(i);
      vec D = normals.col(i);
      int hitIdx = -1; // Initialize to -1 (No Hit)

      for (int j = 0; j < nB; ++j) {
        // Slab AABB test
        double tmin = -INFINITY, tmax = INFINITY;
        bool skip = false;

        for (int k = 0; k < 3; ++k) {
          if (std::abs(D(k)) < 1e-9) {
            if (P(k) < bMins(k, j) || P(k) > bMaxs(k, j)) {
              skip = true; break;
            }
          } else {
            double invD = 1.0 / D(k);
            double t1 = (bMins(k, j) - P(k) -  epsilon) * invD;
            double t2 = (bMaxs(k, j) - P(k) + epsilon) * invD;
            tmin = std::max(tmin, std::min(t1, t2));
            tmax = std::min(tmax, std::max(t1, t2));
          }
        }

        if (skip || tmax < tmin || tmax < 0) continue;

        vec v0 = verticesB.col(facesB(0, j) - 1);
        vec v1 = verticesB.col(facesB(1, j) - 1);
        vec v2 = verticesB.col(facesB(2, j) - 1);

        if (rayIntersectsTriangle(P, D, v0, v1, v2, epsilon) > 0) {
          hitIdx = j + 1;
          break;
        }
      }
      output[i] = hitIdx;
    }
  }
};

// [[Rcpp::export]]
Rcpp::IntegerVector find_first_intersections_parallel(arma::mat centroids, arma::mat normals, arma::mat verticesB, arma::imat facesB, double epsilon) {
  int nA = centroids.n_cols;
  int nB = facesB.n_cols;
  Rcpp::IntegerVector results(nA);

  // Precompute AABBs
  mat bMins(3, nB);
  mat bMaxs(3, nB);
  for (int j = 0; j < nB; ++j) {
    // Armadillo .cols() is faster for extracting the 3 vertices
    uvec idxs = conv_to<uvec>::from(facesB.col(j)) - 1;
    mat tri = verticesB.cols(idxs);
    bMins.col(j) = min(tri, 1);
    bMaxs.col(j) = max(tri, 1);
  }

  IntersectionWorker worker(centroids, normals, verticesB, facesB, bMins, bMaxs, epsilon, results.begin());
  parallelFor(0, nA, worker);
  // worker(0, nA)

  return results;
}


// RcppParallel worker for nearest triangles
struct NearestWorker : public Worker {
  const mat& nodesA;
  const imat& facesA;
  const mat& nodesB;
  const imat& facesB;
  const mat& centroidsB;
  const std::vector<AABB>& boxes;

  // raw pointer to ensure memory is hit
  int* output;

  NearestWorker(const mat& nA, const imat& fA, const mat& nB, const imat& fB, const mat& cB, const std::vector<AABB>& bx, int* res)
    : nodesA(nA), facesA(fA), nodesB(nB), facesB(fB), centroidsB(cB), boxes(bx), output(res) {}

  void operator()(std::size_t begin, std::size_t end) {
    int nA = facesA.n_cols;

    for (std::size_t i = begin; i < end; ++i) {
        uvec idxA = conv_to<uvec>::from(facesA.col(i)) - 1;
        vec3 centroidA = mean(nodesA.cols(idxA), 1);

        double min_dist_squared = datum::inf;
        int best_idx = 0;

        // Iterate through mesh A
        for(int j = 0; j < nA; ++j) {
          // SLAB-inspired pruning: Check if the AABB is even worth looking at
          double d2_box = dist_sq_point_aabb(centroidA, boxes[j]);
          if(d2_box >= min_dist_squared) continue;

          // Actual centroid distance
          vec3 diff = centroidA - centroidsB.col(j);
          double d2 = dot(diff, diff);

          if(d2 < min_dist_squared) {
            min_dist_squared = d2;
            best_idx = boxes[j].face_idx;
          }
        output[i] = best_idx;
      }
    }
  }
};




// [[Rcpp::export]]
Rcpp::IntegerVector find_closest_triangles_parallel(arma::mat nodesA, arma::imat facesA, arma::mat nodesB, arma::imat facesB) {
  int nA = facesA.n_cols;
  int nB = facesB.n_cols;
  IntegerVector results(nA);

  // Precompute mesh B centroids and AABBs for each triangle
  mat centroidsB(3, nB);
  std::vector<AABB> boxes(nB);

  for(int j = 0; j < nB; ++j) {
    uvec idx = conv_to<uvec>::from(facesB.col(j)) - 1;
    mat33 tri = nodesB.cols(idx);

    centroidsB.col(j) = mean(tri, 1);
    boxes[j].min_pt = min(tri, 1);
    boxes[j].max_pt = max(tri, 1);
    boxes[j].face_idx = j + 1;
  }

  NearestWorker worker(nodesA, facesA, nodesB, facesB, centroidsB, boxes, results.begin());
  parallelFor(0, nA, worker);
  // worker(0, nA)

  return results;
}




// [[Rcpp::export]]
Rcpp::IntegerVector find_closest_triangles(arma::mat nodesA, arma::imat facesA, arma::mat nodesB, arma::imat facesB) {
  int nA = facesA.n_cols;
  int nB = facesB.n_cols;
  IntegerVector results(nA);

  // Pre-calculate Mesh B Centroids and AABBs for each triangle
  mat centroidsB(3, nB);
  std::vector<AABB> boxes(nB);

  for(int j = 0; j < nB; ++j) {
    uvec idx = conv_to<uvec>::from(facesB.col(j)) - 1;
    mat33 tri = nodesB.cols(idx);

    centroidsB.col(j) = mean(tri, 1);
    boxes[j].min_pt = min(tri, 1);
    boxes[j].max_pt = max(tri, 1);
    boxes[j].face_idx = j + 1;
  }

  // Iterate through Mesh A
  for(int i = 0; i < nA; ++i) {
    uvec idxA = conv_to<uvec>::from(facesA.col(i)) - 1;
    vec3 centroidA = mean(nodesA.cols(idxA), 1);

    double min_dist_squared = datum::inf;
    int best_idx = 0;

    // Search
    for(int j = 0; j < nB; ++j) {
      // SLAB-inspired pruning: Check if the AABB is even worth looking at
      double d2_box = dist_sq_point_aabb(centroidA, boxes[j]);
      if(d2_box >= min_dist_squared) continue;

      // Actual centroid distance
      vec3 diff = centroidA - centroidsB.col(j);
      double d2 = dot(diff, diff);

      if(d2 < min_dist_squared) {
        min_dist_squared = d2;
        best_idx = boxes[j].face_idx;
      }
    }
    results[i] = best_idx;
  }

  return results;
}



