# polygon.network class----
#' @title S4 class for polygon networks
#' @description Polygon networks are subgraphs made of polygons
#' (i) sharing topographic features and
#' (ii) in contact with the rest of the subgraph by at least 1 polygon edge.
#' Objects of S4 class \code{polygon.network} are typically made using the function
#' \code{\link{poly.network}}
#' @slot membership a vector of numeric values specifying, for each triangle, the
#' index number of the patch to which the triangle is assigned
#' @slot faces a vector of numeric values indicating the mesh triangle indexes
#' @export polygon.network
#' @exportClass polygon.network
polygon.network <- methods::setClass("polygon.network", slots = c(membership = "integer", faces = "integer"))

# poly.network function----
#' @title  Identify polygon networks
#' @description From a selected variable y, identifies patches of adjacent polygons
#' that share a given range of y values. These patches are called ’polygon networks’.
#' @param mesh object of class mesh3d
#' @param y a vector of values to be used to select polygons
#' @param lwr.limit the lower range of values to be selected from y
#' @param upr.limit the upper range of values to be selected from y
#' @param min.size the minimum amount of polygons defining a cluster. Default is set to 3.
#' @return An object of class \code{"\link[=polygon.network-class]{polygon.network}"} composed of
#' the face index and the membership of each triangle answering the set conditions. The function
#' makes patches of contiguous triangles, and each patch is indexed with a unique number corresponding
#' to its membership.
#' @examples
#' #Isolate cusps using elevation:
#' mythreshold <- quantile(elev(dkmodel$cusp), 0.65)
#' cusps <- poly.network(dkmodel$cusp, elev(dkmodel$cusp), lwr.limit = mythreshold,
#' min.size = 100)
#' myvector <- rep(0, Rvcg::nfaces(dkmodel$cusp))
#' myvector[cusps@@faces] <- cusps@@membership[]
#' myvector <- as.factor(myvector)
#' ncusps <- length(levels(myvector)) - 1
#' levels(myvector) <- c(0:ncusps + 1)
#' dkmap(dkmodel$cusp, as.numeric(myvector), col = cbPalette <- c("#000000", "#E69F00",
#' "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
#' col.levels = ncusps + 1, legend.lab = "Elevation (mm)")
#'
#' #Any other variables could be used to define the clusters
#' #Mean curvature:
#' crests <- poly.network(dkmodel$complex, Rvcg::vcgCurve(dkmodel$complex)$meanitmax,
#' lwr.limit = quantile(Rvcg::vcgCurve(dkmodel$complex)$meanitmax, 0.8), min.size = 10)
#' doolkit::dkmap(mesh = dkmodel$complex, y = doolkit::arc(dkmodel$complex,
#' range = c(-20, 20)), col = "arc", col.levels = 256, min.range = -20,
#' max.range = 20, orient = "occlusal", legend.lab = "ARC",
#' alpha.thresh = quantile(doolkit::arc(dkmodel$complex), 0.8), alpha = 0.3,
#' alpha.above = FALSE)
#'
#' valleys <- poly.network(dkmodel$complex, Rvcg::vcgCurve(dkmodel$complex)$meanitmax,
#' upr.limit = quantile(Rvcg::vcgCurve(dkmodel$complex)$meanitmax, 0.2), min.size = 10)
#' doolkit::dkmap(mesh = dkmodel$complex, y = doolkit::arc(dkmodel$complex,
#' range = c(-20, 20)), col = "arc", col.levels = 256, min.range = -20,
#' max.range = 20, orient = "occlusal", legend.lab = "ARC",
#'  alpha.thresh = quantile(doolkit::arc(dkmodel$complex), 0.2), alpha = 0.3,
#'  alpha.above = TRUE)
#'
#' #Orientation and surface of patches:
#' patch_orient <- data.frame(bin = NULL, patch = NULL, size = NULL, surface = NULL)
#' for (i in 1:8) {
#'   Cluster <- poly.network(dkmodel$complex, orient(dkmodel$complex),
#'   lwr.limit = 45 * (i - 1), upr.limit = 45 * i)
#'   Patches <- levels(as.factor(Cluster@@membership))
#'   Bins <- rep(paste(45 * (i - 1), "-", 45 * i), length(Patches))
#'   Areas <- rep(0, length(Patches))
#'   for (j in 1:length(Patches)) {
#'     test <- Cluster@@faces[Cluster@@membership == Patches[j]]
#'     Areas[j] <- round(sum(Rvcg::vcgArea(dkmodel$complex,
#'     perface = TRUE)$pertriangle[test]), 3)
#'   }
#'   patch_orient <- data.frame(rbind(patch_orient,
#'   cbind.data.frame(Bins, Patches, Areas)))
#' }
#' #Since patches made of 3 or less polygons are discarded,
#' #sum of patch areas < total surface area:
#' sum(patch_orient$Areas)
#' Rvcg::vcgArea(dkmodel$complex)
#'
#' @export
poly.network <- function(mesh, y, lwr.limit = stats::quantile(y, 0.75), upr.limit = stats::quantile(y, 1), min.size = 3) {
  # edgelist of faces adjacent to at least another one (excludes borders)
  Edges <- Rvcg::vcgGetEdge(mesh, unique = FALSE)
  Edgelist <- t(matrix(Edges$facept[Edges$border == 0], nrow = 2))
  Edgelist <- data.frame(as.numeric(Edgelist[,1]), as.numeric(Edgelist[,2]))
  # list of faces within the bin
  Facelist = c(1:length(mesh$it[1,]))
  Faces <- Facelist[y >= lwr.limit & y <= upr.limit]
  Faces <- data.frame(Faces)
  # subgraph of faces within the bin AND adjacent to at least another one
  Edgen <- Edgelist[which(Edgelist[, 1] %in% Faces[, 1]), ]
  Edgen <- Edgen[which(Edgen[, 2] %in% Faces[, 1]), ]
  Edgen <- as.matrix(Edgen)
  requireNamespace("igraph", quietly = TRUE)
  Network <- igraph::graph_from_edgelist(Edgen, directed = FALSE)
  # Patches
  Patches <- igraph::components(Network)
  Membership <- Patches$membership[Patches$membership %in% which(Patches$csize >= min.size)]
  # Results
  Results <- methods::new("polygon.network",
                          membership = as.integer(Membership),
                          faces = as.integer(which(Patches$membership %in% Membership)))
  return(Results)
}
