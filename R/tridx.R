# tridx----
#' @title tridx
#' @description Get the pairwise indexing of a 3D mesh relative to another.
#' @param oes object of class mesh3d; should be the outer enamel surface
#' @param edj object of class mesh3d; should be the enamel-dentine junction
#' @param method logical, if TRUE the search is along vertex normals (default is FALSE)
#' @return A data frame.
#' @examples
#' tridx <- data.frame()
#'
#' @export
tridx <- function(oes, edj, method = "closest"){
  # Perform various checks:
  if (!isa(oes, what = "mesh3d")) stop("'oes' must be an object of class 'mesh3d'")
  if (!isa(edj, what = "mesh3d")) stop("'edj' must be an object of class 'mesh3d'")
  if (is.null(oes$normals)) mesh <- Rvcg::vcgUpdateNormals(oes, silent = TRUE)
  if (is.null(edj$normals)) mesh <- Rvcg::vcgUpdateNormals(edj, silent = TRUE)
  if (method != ("closest") || method != ("ortho") || method != ("normal"))

  # Prepare
  nodesA <- oes$vb[1:3, ]
  facesA <- oes$it
  nodesB <- edj$vb[1:3, ]
  facesB <- edj$it
  # Method "closest
  if (method == closest){
    closest_triangles <- find_closest_triangles(nodesA, facesA, nodesB, facesB)
    Results <- data.frame(c(1:length(facesA[1, ])), closest_triangles)
  }

  return(Results)
}
