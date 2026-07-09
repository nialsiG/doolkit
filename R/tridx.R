# tridx----
#' @title tridx
#' @description Get the pairwise indexing of a 3D mesh relative to another.
#' @param oes object of class mesh3d; should be the outer enamel surface
#' @param edj object of class mesh3d; should be the enamel-dentine junction
#' @param method String indicating which method should be used to find triangle pairs.
#' Default is 'nearest', which find the absolute nearest 'edj' triangle.
#' Alternative choices: 'normal' follows the normal to surface A, 'ortho' follows the Z-axis.
#' @param direction Signed integer, indicates the direction for normals / orthogonal methods.
#' Default is '-1' which indicates that mesh B is below mesh A.
#' @param  epsilon Float corresponding to the error margin for intersections.
#' @return A data frame.
#' @examples
#' result <- doolkit::tridx(oes = doolkit::dkmodel$cusp, edj = doolkit::dkmodel$flat, method = "nearest")
#'
#' @export
tridx <- function(oes, edj, method = "nearest", direction = -1, epsilon = 1e-7){
  # Perform various checks:
  if (!isa(oes, what = "mesh3d")) stop("'oes' must be an object of class 'mesh3d'")
  if (!isa(edj, what = "mesh3d")) stop("'edj' must be an object of class 'mesh3d'")
  if (is.null(oes$normals)) mesh <- Rvcg::vcgUpdateNormals(oes, silent = TRUE)
  if (is.null(edj$normals)) mesh <- Rvcg::vcgUpdateNormals(edj, silent = TRUE)
  if (method != ("nearest") && method != ("ortho") && method != ("normal")) stop("'method' must be 'nearest', 'normal' or 'orthogonal'")

  # Prepare
  nodesA <- oes$vb[1:3, ]
  facesA <- oes$it
  nodesB <- edj$vb[1:3, ]
  facesB <- edj$it

  # Computation
  if (method == "nearest"){
    closest_triangles <- find_closest_triangles(nodesA, facesA, nodesB, facesB)
    Results <- data.frame(c(1:length(facesA[1, ])), closest_triangles)
  }
  # nearest_parallel: unstable, causes R session to abort
  else if (method == "nearest_parallel"){
    closest_triangles <- find_closest_triangles_parallel(nodesA, facesA, nodesB, facesB)
    Results <- data.frame(c(1:length(facesA[1, ])), closest_triangles)
  }
  else if (method == "normal"){
    meshA_fixed <- Morpho::facenormals(oes)
    centroids <- meshA_fixed$vb[1:3,]
    normals <- meshA_fixed$normals * direction
    paired_triangles <- find_first_intersections_parallel(centroids, normals, nodesB, facesB, epsilon)
    Results <- data.frame(c(1:length(facesA[1, ])), paired_triangles)
  }
  else if (method == "ortho"){
    centroids <- t(Rvcg::vcgBary(oes))
    ortho_normals <- centroids
    ortho_normals[1, ] <- rep(0, length(centroids[1, ]))
    ortho_normals[2, ] <- rep(0, length(centroids[1, ]))
    ortho_normals[3, ] <- rep(direction, length(centroids[1, ]))
    paired_triangles <- find_first_intersections_parallel(centroids, ortho_normals, nodesB, facesB, epsilon)
    Results <- data.frame(c(1:length(facesA[1, ])), paired_triangles)
  }

  Results[which(Results[, 2] == -1), 2] <- NA

  return(Results)
}
