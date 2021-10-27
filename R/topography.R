# angularity----
#' @title angularity
#' @description Compute the angularity (delta slope).
#' @param mesh object of class mesh3d
#' @param ratio logical, if false standard angularity will be computed (default),
#' if true a relative angularity ratio will be computed (see below)
#' @return If ratio = FALSE, a numeric vector of angularity values i.e. delta slope of
#' each polygon compared to their adjacent polygons, for all the polygons of the
#' mesh.
#' If ratio = TRUE, a numeric vector of angularity ratio values. Ratio is computed
#' by dividing polygon slope by 90, replacing vertex elevation by the average
#' ratio of faces adjacent to the vertex, then dividing the slope of polygons from
#' this new mesh by 90. Although it is a non-standard variable, it was implemented because
#' it better discriminates sharp edges than basic angularity.
#' Warning: both options can take up a significant amount of time for large meshes.
#' @references \href{https://palaeo-electronica.org/2000_1/gorilla/issue1_00.htm}{Ungar and Williamson (2000)}
#' @examples
#' delta_slope <- angularity(dkpongo$OES, ratio = FALSE)
#' summary(delta_slope)
#' #angularity ratio:
#' delta_slope <- angularity(dkpongo$OES, ratio = TRUE)
#' summary(delta_slope)
#'
#' #render on a map:
#' dkmap(dkpongo$OES, angularity(dkpongo$OES), col = "slope", legend.lab = "Angularity (°)")
#' #angularity ratio:
#' dkmap(dkpongo$OES, angularity(dkpongo$OES, ratio = TRUE), col = "slope", legend.lab = "Angularity ratio")
#'
#' @export
angularity <- function(mesh, ratio = FALSE){
  # Perform various checks:
  if (class(mesh) != "mesh3d") stop("mesh must be an object of class 'mesh3d'")
  # Main job
  # ...Get slope (or ratio) of the mesh
  if (ratio == FALSE) Slope <- doolkit::slope(mesh)
  if (ratio == TRUE) Slope <- doolkit::slope(mesh)/90
  # ...replace vertex elevation by the average slope of adjacent faces
  for (i in 1:length(mesh$vb[1,])) mesh$vb[3, i] <- mean(Slope[c(which(mesh$it[1,] == i), which(mesh$it[2,] == i), which(mesh$it[3,] == i))])
  # ...Angularity is the slope (or ratio) of the new mesh
  if (ratio == FALSE) Results <- doolkit::slope(mesh)
  if (ratio == TRUE) Results <- doolkit::slope(mesh)/90
  return(Results)
}

# arc----
#' @title Average-Relative Curvature (ARC)
#' @description Compute a scale-free estimate of mean curvature.
#' @param mesh object of class mesh3d
#' @param range a numeric vector of the form c('minrange', 'maxrange') indicating
#' the set limits for extreme values, beyond which arc values will be truncated.
#' If 'minrange' and 'maxrange' are comprised between 0 and 1, they are used as
#' arc percentages. If 'minrange' and 'maxrange' are not comprised between 0 and 1,
#' they are used as raw arc values
#' @return A numeric vector of area-relative curvature values for all the polygons of the mesh.
#' @details Area-relative curvature (ARC) is obtained by dividing the mean curvature of each triangle
#' by the mean curvature of an hemisphere, the surface area of which is the same as the surface area
#' of the total mesh object. As a result, ARC is a scale-free estimate of surface curvature. It can
#' be used to estimate the sharpness of a mesh.
#' @examples
#' pongo_curvature <- arc(dkpongo$OES)
#' summary(pongo_curvature)
#'
#' #There is a default truncature of extreme values below 1% or above 99%.
#' #Without truncature:
#' pongo_curvature <- arc(dkpongo$OES, range = c(0, 1))
#' summary(pongo_curvature)
#'
#' #mean positive ARC:
#' parc <- mean(pongo_curvature[pongo_curvature >= 0])
#' #mean negative ARC:
#' narc <- mean(pongo_curvature[pongo_curvature < 0])
#'
#' #render on a map:
#' dkmap(dkpongo$OES, arc(dkpongo$OES), col = "arc",
#' min.range = -20, max.range = 20)
#' #absolute truncature of the values above 20 or below -20:
#' dkmap(dkpongo$OES, arc(dkpongo$OES, range = c(-20, 20)), col = "arc",
#' min.range = -20, max.range = 20)
#' @export
arc <- function (mesh, range = c(0.01, 0.99)){
  # Perform various checks:
  if (class(mesh) != "mesh3d") stop("mesh must be an object of class 'mesh3d'")
  # Main job
  # ...Get mean curvature
  Curv <- Rvcg::vcgCurve(mesh)$meanitmax
  # ...Replace border values by 0
  Curv[dkborder(mesh)] <- 0
  # ...Relative truncature
  if (range[1] >= 0 && range[2] <= 1) {
    Curv[Curv > stats::quantile(Curv, range[2])] <- stats::quantile(Curv, range[2])
    Curv[Curv < stats::quantile(Curv, range[1])] <- stats::quantile(Curv, range[1])
  }
  # ...Standardisation
  Curv.theo <- 2.4481 * (1 / sqrt(Rvcg::vcgArea(mesh)))
  Arc <- Curv/Curv.theo
  # ...Absolute truncature
  if (range[1] < 0) Arc[Arc < range[1]] <- range[1]
  if (range[2] > 1) Arc[Arc > range[2]] <- range[2]
  return(Arc)
}

# area2d----
#' @title 2D surface area
#' @description Compute the area of a 2d projection on the (xy) plane.
#' @param mesh object of class mesh3d
#' @param method the method used to compute the hull of the 2d projection, either
#' 'convex' (see \code{\link[grDevices]{chull}}),
#' 'concave' (see \code{\link[concaveman]{concaveman}})
#' or 'alpha' (see \code{\link[alphahull]{ahull}}). The default method is 'concave'.
#' @param alpha the value of alpha, between 0 and 1; see \code{\link[alphahull]{ahull}}.
#' Only used if method = 'alpha'.
#' @return A single 2D surface area value, corresponding to the footprint of the mesh.
#' @seealso \code{\link{rfi}}
#' @examples
#' #The following objects should have the exact same footprints:
#' area2d(dkmodel$basin)
#' area2d(dkmodel$complex)
#' area2d(dkmodel$cusp)
#' area2d(dkmodel$flat)
#'
#' #Graphical rendering of convex hull
#' x <- dkpongo$OES
#' FootprintVerts <- t(x$vb[1:2, ])
#' Hull <- grDevices::chull(x = FootprintVerts[, 1], y = FootprintVerts[, 2])
#' plot(x$vb[1, ], x$vb[2, ], xlab = "x", ylab = "y")
#' points(x$vb[1, Hull], x$vb[2, Hull], col = "orange1")
#'
#' #Graphical rendering of concave hull
#' x <- dkpongo$OES
#' FootprintVerts <- t(x$vb[1:2, ])
#' FootprintVerts[, 1] <- FootprintVerts[, 1] - min(FootprintVerts[, 1])
#' FootprintVerts[, 2] <- FootprintVerts[, 2] - min(FootprintVerts[, 2])
#' Hull <- concaveman::concaveman(points = FootprintVerts)
#' plot(x$vb[1, ] - min(x$vb[1, ]), x$vb[2, ] - min(x$vb[2, ]), xlab = "x", ylab = "y")
#' points(Hull, col = "green1")
#'
#' #Graphical rendering of alpha hull for alpha = 0.1
#' x <- dkpongo$OES
#' FootprintVerts <- t(x$vb[1:2, ])
#' Hull <- alphahull::ahull(x = FootprintVerts[, 1], y = FootprintVerts[, 2], alpha = 0.1)
#' plot(x$vb[1, ], x$vb[2, ], xlab = "x", ylab = "y")
#' points(x$vb[1, Hull$ashape.obj$edges[, 1:2]], x$vb[2, Hull$ashape.obj$edges[, 1:2]], col = "cyan1")
#' @export
area2d <- function(mesh, method = "concave", alpha = 0.5){
  # Perform various checks:
  if (class(mesh) != "mesh3d") stop("mesh must be an object of class 'mesh3d'")
  # Main job
  # ...Get the (x,y) coordinates of all the vertices
  FootprintVerts <- t(mesh$vb[1:2, ])
  FootprintVerts <- FootprintVerts[!duplicated(FootprintVerts),]
  if (method == "alpha") {
    #1 extract the hull using the package aphahull
    Hull <- alphahull::ahull(x = FootprintVerts[, 1], y = FootprintVerts[, 2], alpha = alpha)
    #2 get the 2D area using the package aphahull
    Area2D <- alphahull::areaahull(Hull)
  }
  if (method == "concave") {
    #1 move all coordinates >= 0
    FootprintVerts[, 1] <- FootprintVerts[, 1] - min(FootprintVerts[, 1])
    FootprintVerts[, 2] <- FootprintVerts[, 2] - min(FootprintVerts[, 2])
    #2 extract the hull using the package aphahull
    Hull <- concaveman::concaveman(points = FootprintVerts)
    #3 make a polygon from the points of the hull and measure its area
    Poly <- sp::Polygon(Hull, hole = FALSE)
    Area2D <- Poly@area
  }
  if (method == "convex") {
    #1 extract the hull using chull from grDevices
    Hull <- grDevices::chull(x = FootprintVerts[, 1], y = FootprintVerts[, 2])
    #2 close the contour
    FootprintVerts <- rbind(FootprintVerts, FootprintVerts[1, ])
    #3 make a polygon from the points of the hull and measure its area
    Poly <- sp::Polygon(FootprintVerts[Hull, ], hole = FALSE)
    Area2D <- Poly@area
  }
  return(Area2D)
}

# dne----
#' @title Dirichlet normal energy
#' @description Compute the Dirichlet normal energy.
#' @param mesh object of class mesh3d
#' @param range an integer between 0 and 1 indicating the percentage of values
#' to consider for the computation. Following Pampush et al. (2016) default is
#' set to 0.999.
#' @param total logical, should the result of the function be the total DNE.
#' @return If total = FALSE, a numeric vector of dne values for all the polygons of the
#' mesh. If total = TRUE, a single DNE value, calculated as the sum of the
#' products of polygon normal energies * polygon areas.
#' @references \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/ajpa.21489}{Bunn et al. (2011)}
#' @references \href{https://link.springer.com/article/10.1007/s10914-016-9326-0}{Pampush et al. (2016)}
#' @seealso \code{\link[molaR]{molaR::DNE}}
#' @examples
#' pongo_dne <- dne(dkpongo$OES)
#' summary(pongo_dne)
#'
#' #total DNE value corresponds to the sum of products Dne * Area:
#' round(sum(dne(dkpongo$OES)*Rvcg::vcgArea(dkpongo$OES, perface = TRUE)$pertriangle), 3)
#' #can be directly computed using \code{dne}:
#' dne(dkpongo$OES, total = TRUE)
#'
#' #render on a map:
#' dkmap(dkpongo$OES, dne(dkpongo$OES), legend.type = "log", col = "dne")
#' @export
dne <- function(mesh, range = 0.999, total = FALSE){
  # Perform various checks:
  if (class(mesh) != "mesh3d") stop("mesh must be an object of class 'mesh3d'")
  if (is.null(mesh$normals)) mesh <- Rvcg::vcgUpdateNormals(mesh, silent = TRUE)
  # Main job
  Dne <- NULL
  for (i in 1:length(mesh$it[1, ])){
    N1T <- mesh$vb[1:3, mesh$it[1, i]]
    N2T <- mesh$vb[1:3, mesh$it[2, i]]
    N3T <- mesh$vb[1:3, mesh$it[3, i]]
    NN1T <- mesh$normals[, mesh$it[1, i]]
    NN2T <- mesh$normals[, mesh$it[2, i]]
    NN3T <- mesh$normals[, mesh$it[3, i]]
    U <- c(N2T[1]-N1T[1], N2T[2]-N1T[2], N2T[3]-N1T[3])
    V <- c(N3T[1]-N1T[1], N3T[2]-N1T[2], N3T[3]-N1T[3])
    NU <- c(NN2T[1]-NN1T[1], NN2T[2]-NN1T[2], NN2T[3]-NN1T[3])
    NV <- c(NN3T[1]-NN1T[1], NN3T[2]-NN1T[2], NN3T[3]-NN1T[3])
    UU <- U%*%U
    UV <- U%*%V
    VV <- V%*%V
    NUNU <- NU%*%NU
    NUNV <- NU%*%NV
    NVNV <- NV%*%NV
    G <- matrix(c(UU, UV, UV, VV), nrow = 2, ncol = 2, byrow = TRUE)
    H <- matrix(c(NUNU, NUNV, NUNV, NVNV), nrow = 2, ncol = 2, byrow = TRUE)
    Dne[i] <- sum(diag(MASS::ginv(G)%*%H))
  }
  # ...Set border values to 0
  Dne[doolkit::dkborder(mesh)] <- 0
  # ...Set values out of range to 0
  Dne[Dne >= stats::quantile(Dne, range)] <- 0
  # ...Return vector (if total = FALSE) or single value (if total = TRUE)
  if(!total) return(Dne)
  if(total) return(round(sum(Dne * Rvcg::vcgArea(mesh, perface = TRUE)$pertriangle), 3))
}

# elev----
#' @title elevation
#' @description Compute the elevation (z component of triangle barycenter).
#' @param mesh object of class mesh3d
#' @param origin logical, if TRUE the z of the mesh is adjusted so that the lowest z = 0
#' (see \code{\link{dkorigin}})
#' @return A numeric vector of elevation values for all the polygons of the mesh.
#' @seealso \code{\link{inclin}}
#' @seealso \code{\link{rfi}}
#' @seealso \code{\link{slope}}
#' @examples
#' pongo_elev <- elev(dkpongo$OES)
#' summary(pongo_elev)
#'
#' #render on a map:
#' dkmap(dkpongo$OES, elev(dkpongo$OES))
#' @export
elev <- function(mesh, origin = TRUE){
  # Perform various checks:
  if (class(mesh) != "mesh3d") stop("mesh must be an object of class 'mesh3d'")
  # Main job
  # ...Set back origin to 0 if origin is TRUE:
  if (origin) {
    mesh <- doolkit::dkorigin(mesh)
  }
  # ...computation
  Elev <- Rvcg::vcgBary(mesh)[, 3]
  return(Elev)
}

# hypso----
#' @title hypso
#' @description Compute the maximum height, length, width and corresponding
#' hypsodonty index (ratio of the maximum height over the maximum length)
#' @param mesh object of class mesh3d
#' @param origin logical, whether to set the z of the mesh's lowermost point
#' to zero.
#' @return A list of values for hypsodonty index, height, length and width of
#' the mesh. The hypsodonty index is not expressed relative to 100 but to 1.
#'     Note: the tooth surface is expected to be oriented such as the X-axis
#'     is the bucco-lingual axis, the Y-axis is the mesio-distal axis, and the
#'     occlusal plane is parallel to the (XY) plane and faces upward.
#' @examples
#' hypso(dkpongo$OES)
#' @export
hypso <- function(mesh, origin = TRUE){
  # Perform various checks:
  if (class(mesh) != "mesh3d") stop("mesh must be an object of class 'mesh3d'")

  # Main job
  # ...Set back origin to 0
  if (origin) {
    mesh <- doolkit::dkorigin(mesh)
  }
  # ...Get x, y and z dimensions
  Xmax <- max(mesh$vb[1,])
  Xmin <- min(mesh$vb[1,])
  Ymax <- max(mesh$vb[2,])
  Ymin <- min(mesh$vb[2,])
  Zmax <- max(mesh$vb[3,])
  Zmin <- min(mesh$vb[3,])
  # ...Computation
  Hypso <- list()
  Hypso$hypsodonty <- round((Zmax - Zmin) / (Ymax - Ymin), 2)
  Hypso$height <- round(Zmax - Zmin, 2)
  Hypso$length <- round(Ymax - Ymin, 2)
  Hypso$width <- round(Xmax - Xmin, 2)
  return(Hypso)
}

# inclin----
#' @title inclin
#' @description Compute inclination i.e. the angle between triangles and the vertical plane in degrees, comprised
#' between 0 and 180.
#' @param mesh object of class mesh3d
#' @return A numeric vector of inclination values for all the polygons of the mesh.
#' @references \href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0066142}{Guy et al. (2013)}
#' @seealso \code{\link{slope}}
#' @examples
#' pongo_inclin <- inclin(dkpongo$OES)
#' summary(pongo_inclin)
#'
#' #render on a map:
#' dkmap(dkpongo$OES, inclin(dkpongo$OES), col = "inclin",
#' min.range = 0, max.range = 180, legend = TRUE)
#' @export
inclin <- function(mesh){
  # Perform various checks:
  if (class(mesh) != "mesh3d") stop("mesh must be an object of class 'mesh3d'")

  # Main job
  # ...Get normals
  Normals <- Morpho::facenormals(mesh)$normals
  # ...Calculate inclination from the angle of normals
  NU <- sqrt((Normals[1, ])^2 + (Normals[2, ])^2)
  NV <- sqrt((Normals[1, ])^2 + (Normals[2, ])^2 + (Normals[3, ])^2)
  Rad <- acos(NU/NV)
  Deg <- 180*(Rad)/pi
  Inclin <- Deg
  # ...Correct for polygons with negative normals = facing downward
  Inclin[Normals[3, ] < 0] <- 90-Deg[Normals[3, ] < 0]
  Inclin[Normals[3, ] >= 0] <- 90+Deg[Normals[3, ] >= 0]
  return(Inclin)
}

# oedist----
#' @title Distance from outer enamel surface to enamel dentine junction
#' @description Compute the distance from enamel vertices to dentine mesh.
#' @param oes object of class mesh3d; should be the outer enamel surface
#' @param edj object of class mesh3d; should be the enamel-dentine junction
#' @param ray logical, if TRUE the search is along vertex normals (default)
#' @return A numeric vector of vertex-to-mesh distance values for all the polygons of the x mesh.
#' @references \href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0066142}{Guy
#' et al. (2013)}
#' @references \href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0138802}{Guy
#' et al. (2015)}
#' @references \href{https://www.frontiersin.org/articles/10.3389/fphys.2017.00524/full}{Thiery
#' et al. (2017)}
#' @references \href{https://royalsocietypublishing.org/doi/10.1098/rsbl.2019.0671}{Schwartz
#' et al. (2020)}
#' @seealso \code{\link[Morpho]{Morpho::meshDist}}
#' @examples
#' edd_pongo <- oedist(dkpongo$OES, dkpongo$EDJ)
#' summary(edd_pongo)
#' AETgeom <- mean(edd_pongo)
#' #Geometric relative enamel thickness, obtained by dividing AETgeom by the square root of EDJ area
#' #Note: it is different from classic RET which requires the volume of the dentine inside the enamel cap (see Thiery et al., 2017)
#' AETgeom/sqrt(Rvcg::vcgArea(dkpongo$EDJ))
#' #Absolute crown strength:
#' edj_radius_pongo <- max(dist(cbind(dkpongo$EDJ$vb[1,], dkpongo$EDJ$vb[2,])))/2
#' sqrt(mean(edd_pongo) * edj_radius_pongo)
#'
#' #render on a map:
#' dkmap(dkpongo$OES, oedist(dkpongo$OES, dkpongo$EDJ))
#' #distance map can also be rendered on EDJ surface:
#' dkmap(dkpongo$EDJ, oedist(dkpongo$EDJ, dkpongo$OES, ray = FALSE))
#' @export
oedist <- function(oes, edj, ray = TRUE){
  # Perform various checks:
  if (class(oes) != "mesh3d") stop("oes must be an object of class 'mesh3d'")
  if (class(edj) != "mesh3d") stop("edj must be an object of class 'mesh3d'")
  # Main job
  # ...Compute absolute shortest distances FROM vertex TO second mesh
  VertDist <- abs(Morpho::meshDist(oes, edj, distvec = NULL, sign = FALSE, plot = FALSE, ray = ray)$dists)
  # ...Get face value = average of corresponding vertex values
  TrDist <- rep(0, length(oes$it[1, ]))
  for (i in 1:length(oes$it[1, ])) {
    Dist1 <- VertDist[oes$it[1, i]]
    Dist2 <- VertDist[oes$it[2, i]]
    Dist3 <- VertDist[oes$it[3, i]]
    TrDist[i] <- mean(c(Dist1, Dist2, Dist3))
  }
  return(TrDist)
}

# opc----
#' @title orientation patch count
#' @description Count the number of orientation patches using \code{\link[poly.network]{poly.network}}.
#' @param mesh object of class mesh3d
#' @param bins the number of orientation bins to be defined (default set to 8)
#' @param min.size the minimal amount of polygons defining a "patch" (default set to 3)
#' @param rotation if applicable, the number of degrees to which bins are to be rotated.
#' By default the bins start from an angle of pi/2 and rotates clockwise.
#' @return A data.frame displaying the number of patches and their size (number of triangles)
#' for each orientation bin. Note: if you want the surface area of each patch, see\code{\link[poly.network]{poly.network}}
#' @references \href{https://www.nature.com/articles/nature05433}{Evans et al. (2007)}
#' @seealso \code{\link{orient}}
#' @seealso \code{\link{opcr}}
#' @examples
#' #8 bins (default):
#' pongo_opc <- opc(dkpongo$OES)
#' #8 bins starting from mesial, as in Evans et al. 2007:
#' pongo_opc <- opc(dkpongo$OES, rotation = -(360/16))
#' #4 bins (mesial, buccal, distal and lingual):
#' pongo_opc <- opc(dkpongo$OES, bins = 4, rotation = -(360/8))
#'
#' @export
opc <- function(mesh, bins = 8, min.size = 3, rotation = 0){
  # Perform various checks:
  if (class(mesh) != "mesh3d") stop("mesh must be an object of class 'mesh3d'")
  # Main job
  # ...A function to count patches:
  patchcount <- function(mesh, bins, min.size, rotation) {
    #1: get orientation
    d <- data.frame(face = c(1:length(mesh$it[1,])), orient = doolkit::orient(mesh))
    #2: get edgelist of faces adjacent to at least another one (excludes borders)
    Edges <- Rvcg::vcgGetEdge(mesh, unique = FALSE)
    Edgelist <- t(matrix(Edges$facept[Edges$border == 0], nrow = 2))
    Edgelist <- data.frame(as.numeric(Edgelist[,1]), as.numeric(Edgelist[,2]))
    #3: get patch-related data for each bin
    Patch.data <- data.frame(Bin = 0, Count = 0, Size = 0)
    for (i in (1:bins)) {
      bin_low_limit <- (i-1)*(360/bins) + rotation
      bin_up_limit <- i*(360/bins) + rotation
      # get list of faces within the bin
      Faces <- d$face[d$orient >= bin_low_limit & d$orient < bin_up_limit]
      if(bin_up_limit >= 360) {
        bin_up_limit <- bin_up_limit - 360
        Faces2 <- d$face[d$orient < bin_up_limit]
        Faces <- data.frame(c(Faces, Faces2))
        if (bin_low_limit == 360) bin_low_limit <- 0
      }
      Faces <- data.frame(Faces)
      # for a given bin, get subgraph of the faces adjacent to at least another one
      edgen <- Edgelist[which(Edgelist[, 1] %in% Faces[, 1]), ]
      edgen <- edgen[which(edgen[, 2] %in% Faces[, 1]), ]
      edgen <- as.matrix(edgen)
      requireNamespace("igraph", quietly = TRUE)
      Network <- igraph::graph_from_edgelist(edgen, directed = FALSE)
      # get clusters within the bin
      Patches <- igraph::components(Network)
      Patch.data[i, 1] <- paste(bin_low_limit, "-", bin_up_limit, sep = "")
      Patch.data[i, 2] <- length(Patches$csize[Patches$csize >= min.size])
      Patch.data[i, 3] <- sum(Patches$csize[Patches$csize >= min.size])
    }
    return(Patch.data)
  }
  # ...Use function to get patches
  Opc <- list()
  Opc$patches <- patchcount(mesh = mesh, bins = bins, min.size = min.size, rotation = rotation)
  Opc$opc <- sum(Opc$patches$Count)
  return(Opc)
}

# opcr----
#' @title orientation patch count rotated
#' @description Compute the orientation patch count rotated of a triangle mesh.
#' @param mesh object of class mesh3d
#' @param bins the number of orientation bins to be defined (default set to 8)
#' @param min.size the minimal amount of polygons defining a "patch" (default set to 3)
#' @return A data.frame displaying the number of patches and their size (number of triangles)
#' for each orientation bin.
#' @references \href{https://www.nature.com/articles/nature10880}{Wilson et al. (2012)}
#' @seealso \code{\link{opc}}
#' @seealso \code{\link{orient}}
#' @seealso \code{\link[molaR]{OPCr}}
#' @examples
#' #8bins (default):
#' pongo_opcr <- opcr(dkpongo$OES)
#' #16 bins (computation time increase exponentially):
#' pongo_opcr <- opcr(dkpongo$OES, bins = 16)
#' @export
opcr <- function(mesh, bins = 8, min.size = 3){
  # Perform various checks:
  if (class(mesh) != "mesh3d") stop("mesh must be an object of class 'mesh3d'")
  # Main job
  Opcr <- list(opcr = 0, rotations = data.frame(Rotation = 0, Opc = 0))
  for (j in (1:bins)) {
    Rotated.angle <- (j-1)*(360/(bins*bins))
    Patches <- doolkit::opc(mesh = mesh, bins = bins, min.size = min.size, rotation = Rotated.angle)$patches
    Opcr$rotations[j, 1] <- paste(Rotated.angle)
    Opcr$rotations[j, 2] <- sum(Patches$Count)
  }
  Opcr$opcr <- mean(Opcr$rotations$Opc)
  return(Opcr)
}

# orient----
#' @title orientation of polygons
#' @description Returns the occlusal orientation (exposure in GIS)
#' @param mesh object of class mesh3d
#' @return A numeric vector of occlusal orientation values in degrees for all
#' the polygons of the mesh.
#' Let the orientation from above be depicted as a trigonomical circle,
#' then for a tooth positioned as in Guy et al. (2015) an orientation of 0 (mesial) would
#' be located at an angle of pi/2, and an orientation of 90° (buccal) would
#' be located at an angle of 2*pi.
#' @seealso \code{\link{opc}}
#' @seealso \code{\link{opcr}}
#' @examples
#' pongo_orient <- orient(dkpongo$OES)
#' summary(pongo_orient)
#'
#' #render on a map:
#' #2 bins
#' dkmap(dkpongo$OES, orient(dkpongo$OES), col.levels = 2, col = c("yellow1", "red1"),
#' min.range = 0, max.range = 360, legend = TRUE, legend.type = "pie")
#' #4 bins
#' dkmap(dkpongo$OES, orient(dkpongo$OES), col.levels = 4, col = c("royalblue1",
#' "yellow1", "red1", "green3"), min.range = 0, max.range = 360, legend = TRUE,
#' legend.type = "pie")
#' #8 bins
#' dkmap(dkpongo$OES, orient(dkpongo$OES), col.levels = 8, col = "orient",
#' min.range = 0, max.range = 360, legend = TRUE, legend.type = "pie")
#' @export
orient <- function(mesh){
  # Perform various checks:
  if (class(mesh) != "mesh3d") stop("mesh must be an object of class 'mesh3d'")
  # Main job
  # ...Get face normals
  Normals <- Morpho::facenormals(mesh)$normals
  # ...Get face orientation
  memb1 <- atan2(Normals[1,], Normals[2,])
  memb1[memb1<0] <- memb1[memb1<0]+(2*pi)
  memb2 <- atan2(1,0)
  theta <- 180*(memb1 - memb2)/pi+90
  return(theta)
}

# rfi----
#' @title relief index
#' @description Compute the relief index of a  3d triangle mesh.
#' @param mesh object of class mesh3d
#' @param method a character string indicating which method is to be used for the computation of
#' relief index
#' @param hull the method used to compute the hull of the 2d projection, either 'convex' 'concave' or 'alpha'.
#' The default method is 'convex'. See \code{\link{area2d}}
#' @param alpha the alpha value when using alpha hull to compute 2d area. See \code{\link[alphahull]{ahull}}
#' @return A single relief index value.
#' @references \href{https://doi.org/10.1016/j.jhevol.2008.08.002}{Boyer (2008)}
#' \href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0066142}{Guy
#' et al. (2013)}
#' \href{https://palaeo-electronica.org/2000_1/gorilla/issue1_00.htm}{Ungar and Williamson (2000)}
#' @seealso \code{\link{area2d}}
#' @seealso \code{\link[molaR]{RFI}}
#' @examples
#' pongo_rfi <- rfi(dkpongo$OES, method = "Ungar")
#' pongo_lrfi <- rfi(dkpongo$OES, method = "Boyer")
#' pongo_gamma <- rfi(dkpongo$OES, method = "Guy")
#' @export
rfi <- function(mesh, method = "Ungar", hull = "concave", alpha = 0.5){
  # Perform various checks:
  if (class(mesh) != "mesh3d") stop("mesh must be an object of class 'mesh3d'")
  # Main job
  if (method == "Ungar"){
    Surf3D <- Rvcg::vcgArea(mesh)
    Surf2D <- doolkit::area2d(mesh, method = hull, alpha = alpha)
    RFI <- Surf3D/Surf2D
  }
  if (method == "Boyer"){
    Surf3D <- Rvcg::vcgArea(mesh)
    Surf2D <- doolkit::area2d(mesh, method = hull, alpha = alpha)
    RFI <- log(sqrt(Surf3D)/sqrt(Surf2D))
  }
  if (method == "Guy"){
    Dataset <- data.frame(Area = Rvcg::vcgArea(mesh), Elevation = doolkit::elev(mesh), Slope = doolkit::slope(mesh))
    RFI <- (sum(Dataset$Area[Dataset$Slope<45])*mean(Dataset$Elevation[Dataset$Slope<45])) / (sum(Dataset$Area[Dataset$Slope>45])*mean(Dataset$Elevation[Dataset$Slope>45]))
  }
  return(RFI)
}

# rrate----
#' @title relief rate
#' @description Compute the relief rate from a sub-sample of a 3d triangle mesh.
#' For instance, the relief rate could be computed from the portion of a molar above
#' the lowermost point of its central basin, compared to the whole tooth.
#' @param uncropped object of class mesh3d. Should entirely contain the 'cropped' argument.
#' @param cropped object of class mesh3d. Should be part of the 'uncropped' argument.
#' @param origin logical, if TRUE both cropped and uncropped mesh are translated along the z-axis
#' so that the lowest z of the uncropped mesh = 0; see \code{\link{dkorigin}}
#' @return A single relief rate value.
#' @references \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/ajpa.23916}{Thiery et al. (2019)}
#' @examples
#' pongo_medelev <- median(elev(dkpongo$OES))
#' pongo_basins <- dkcrop(dkpongo$OES, which(elev(dkpongo$OES) < pongo_medelev))
#' pongo_cusps <- dkcrop(dkpongo$OES, which(elev(dkpongo$OES) > pongo_medelev))
#'
#' #Not run:
#' dkmap(pongo_basins, elev(pongo_basins))
#' dkmap(pongo_cusps, elev(pongo_cusps))
#'
#' rrate(dkpongo$OES, pongo_basins)
#' rrate(dkpongo$OES, pongo_cusps)
#' @export
rrate <- function(uncropped, cropped, origin = TRUE){
  # Perform various checks:
  if (class(uncropped) != "mesh3d") stop("uncropped must be an object of class 'mesh3d'")
  if (class(cropped) != "mesh3d") stop("cropped must be an object of class 'mesh3d'")
  # Main job
  # ...Set back origin to 0
  if (origin) {
    cropped$vb[3, ] <- cropped$vb[3, ] - min(uncropped$vb[3, ])
    uncropped$vb[3, ] <- uncropped$vb[3, ] - min(uncropped$vb[3, ])
  }
  # ...computation
  RRate <- mean(elev(cropped) - min(elev(cropped))) / mean(elev(uncropped) - min(elev(uncropped)))
  return (RRate)
}

# shape.index----
#' @title shape.index
#' @description Compute various shape indices.
#' @param mesh object of class mesh3d
#' @param origin logical, if TRUE the z of the mesh is adjusted so that the lowest z = 0;
#' see \code{\link{dkorigin}}
#' @return A list of indices:
#' \itemize{
#'   \item First item
#'   \item Second item
#'   \item Second item
#'   \item Second item
#' }
#' @details A handful of indices have been developed to characterize the shape of natural
#' landscapes, including drainage basins. While some indices are very scale-sensitive (e.g.,
#' Gravelius' compactness coefficient, see Bendjoudi & Hubert, 2002), others are dimensionless.
#' Horton (1932) introduced a form factor computed as the quotient of the basin's surface area
#' over the square of the maximum basin length. Schumm (1956) developed a basin elongation index
#' computed as the quotient of twice the square root of surface area over the product of basin
#' length and the squareroot of pi. Lastly, Chorley et al. (1957) developed a lemniscate ratio
#' which corresponds to the ratio between the surface of a lemniscate of same length over the
#' basin area,and computed as (pi*(Length^2))/(4*Area).
#' @references \href{https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/TR013i001p00350}{Horton (1932)}
#' @references \href{https://pubs.geoscienceworld.org/gsa/gsabulletin/article-abstract/67/5/597/4811/EVOLUTION-OF-DRAINAGE-SYSTEMS-AND-SLOPES-IN?redirectedFrom=fulltext}{Schumm (1956)}
#' @references \href{http://www.ajsonline.org/content/255/2/138.extract}{Chorley et al. (1957)}
#' @references \href{http://hydrologie.org/hsj/470/hysj_47_06_0921.pdf}{Bendjoudi & Hubert (2002)}
#' @examples
#' shape.index(dkpongo$OES)
#' @export
shape.index <- function(mesh, origin = TRUE) {
  # Perform various checks:
  if (class(mesh) != "mesh3d") stop("mesh must be an object of class 'mesh3d'")
  # Main job
  # ...Set back origin to 0
  if (origin) {
    mesh <- doolkit::dkorigin(mesh)
  }
  # ...Get area and length
  Area <- Rvcg::vcgArea(mesh)
  Length <- hypso(mesh)$length
  # ...Calculate indices
  Indices <- list(FormFactor = NULL, Elongation = NULL, K = NULL)
  Indices$FormFactor <- Area/(Length^2)
  Indices$Elongation <- (2*sqrt(Area))/(Length*sqrt(pi))
  Indices$K <- (pi*(Length^2))/(4*Area)
  return(Indices)
}

# slope----
#' @title slope
#' @description Compute slope i.e. the angle between triangles and the horizontal plane in degrees, comprised
#' between 0 and 90.
#' @param mesh object of class mesh3d
#' @return A numeric vector of slope values for all the polygons of the mesh.
#' @seealso \code{\link{inclin}}
#' @examples
#' pongo_slope <- slope(dkpongo$OES)
#' summary(pongo_slope)
#'
#' #render on a map:
#' dkmap(dkpongo$OES, slope(dkpongo$OES), col.levels = 9, col = "slope",
#' min.range = 0, max.range = 90, legend = TRUE)
#' @export
slope <- function(mesh) {
  # Perform various checks:
  if (class(mesh) != "mesh3d") stop("mesh must be an object of class 'mesh3d'")
  # Main job
  # ...Get surface normals
  Normals <- Morpho::facenormals(mesh)$normals
  # ...Calculate slope from the angle of normals
  NU <- sqrt((Normals[1, ])^2 + (Normals[2, ])^2)
  NV <- sqrt((Normals[1, ])^2 + (Normals[2, ])^2 + (Normals[3, ])^2)
  Rad <- acos(NU/NV)
  Deg <- 180*(Rad)/pi
  Slope <- 90-Deg
  return(Slope)
}
