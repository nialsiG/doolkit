#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Helper function to avoid calling the BLAS 'dot'
double manual_dot(const vec& a, const vec& b) {
  double res = 0;
  for(uword i = 0; i < a.n_elem; ++i) {
    res += a[i] * b[i];
  }
  return res;
}

// [[Rcpp::export]]
IntegerVector find_closest_triangles(arma::mat nodesA, arma::imat facesA, arma::mat nodesB, arma::imat facesB) {
  int nA = facesA.n_cols;
  int nB = facesB.n_cols;
  IntegerVector closest_indices(nA);

  // Pre-calculate centroids for Mesh B
  arma::mat centroidsB(3, nB);
  for(int j = 0; j < nB; ++j) {
    // R to C++ indexing: subtract 1 from the face indices
    uvec idx = conv_to<uvec>::from(facesB.col(j)) - 1;
    centroidsB.col(j) = (nodesB.col(idx(0)) + nodesB.col(idx(1)) + nodesB.col(idx(2))) / 3.0;
  }
  // Iterate through Mesh A triangles
  for(int i = 0; i < nA; ++i) {
    uvec idxA = conv_to<uvec>::from(facesA.col(i)) - 1;
    vec centroidA = (nodesA.col(idxA(0)) + nodesA.col(idxA(1)) + nodesA.col(idxA(2))) / 3.0;

    double min_dist_sq = datum::inf; // Use squared distance for speed
    int best_idx = 0;

    // Find closest centroid in B
    for(int j = 0; j < nB; ++j) {
      vec diff = centroidA - centroidsB.col(j);
      double d2 = manual_dot(diff, diff);

      if(d2 < min_dist_sq) {
        min_dist_sq = d2;
        best_idx = j + 1; // Return 1-based index for R
      }
    }
    closest_indices[i] = best_idx;
  }

  return closest_indices;
}



