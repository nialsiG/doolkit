# dkborder----
#' @title dkborder
#' @description Selects the border of a 3d triangle mesh.
#' @param mesh object of class mesh3d
#' @return A vector of indices corresponding to the triangles with at least one vertex on the border of the mesh.
#' @examples
#' border <- dkborder(dkmodel$cusp)
#'
#' # Map the border in orange:
#' is_border <- rep(1, Rvcg::nfaces(dkmodel$cusp))
#' is_border[border] <- 2
#' dkmap(dkmodel$cusp, is_border, col = c("white", "#E69F00"), col.levels = 2, legend = FALSE,
#' scalebar = FALSE, smooth = FALSE)
#'
#' # Compare with vcgBorder from the R package Rvcg, in blue:
#' vcgborder <- which(Rvcg::vcgBorder(dkmodel$cusp)$borderit == TRUE)
#' is_border[vcgborder] <- 3
#' dkmap(dkmodel$cusp, is_border, col = c("white", "#E69F00", "#56B4E9"), col.levels = 3,
#' legend = FALSE, scalebar = FALSE, smooth = FALSE)
#' #As you can see, it all depends on what you want to select!
#'
#' @export
dkborder <- function(mesh){
  # Perform various checks:
  if (!isa(mesh, what = "mesh3d")) stop("mesh must be an object of class 'mesh3d'")
  # Main job
  # ...Make data.frame with face number and triangle vertices
  Mesh <- data.frame(face = c(1:length(mesh$it[1,])), t(mesh$it))
  # ...Get the faces with vertices on the border
  Bordervb <- which(Rvcg::vcgBorder(mesh)$bordervb == TRUE)
  Match1 <- Mesh$face[which(Mesh[, 2] %in% Bordervb)]
  Match2 <- Mesh$face[which(Mesh[, 3] %in% Bordervb)]
  Match3 <- Mesh$face[which(Mesh[, 4] %in% Bordervb)]
  # ...Return the border
  TrueBorderit <- c(Match1, Match2, Match3)
  return(TrueBorderit)
}

# dkcrop----
#' @title crop a mesh
#' @description Crop a 3d triangle mesh.
#' @param mesh object of class mesh3d
#' @param y numeric vector indicating which polygons should be cropped; or an object of class \code{polygon.network}
#' @return A new object of class mesh3d for which all polygons out of y have been removed.
#' @examples
#' #Crop above a certain threshold:
#' mythreshold <- quantile(elev(dkmodel$basin), 0.5)
#' mypolynetwork <- poly.network(dkmodel$basin, elev(dkmodel$basin),
#' lwr.limit = mythreshold)
#' mynewmesh <- dkcrop(dkmodel$basin, mypolynetwork)
#' dkmap(mynewmesh, elev(mynewmesh))
#'
#' #Crop the sharpest dental elements:
#' sharpmesh <- dkcrop(dkmodel$basin, poly.network(dkmodel$basin,
#' Rvcg::vcgCurve(dkmodel$basin)$meanitmax,
#' lwr.limit = quantile(Rvcg::vcgCurve(dkmodel$basin)$meanitmax, 0.8),
#' min.size = 50))
#' dkmap(sharpmesh, arc(sharpmesh), col = "arc", col.levels = 20,
#' min.range = -20, max.range = 20)
#' #Map of the sharpest elements' elevation, slope and orientation;
#' dkmap(sharpmesh, elev(sharpmesh), col = "elev")
#' dkmap(sharpmesh, slope(sharpmesh), col = "slope", col.levels = 9,
#' min.range = 0, max.range = 90)
#' dkmap(sharpmesh, orient(sharpmesh), col = "orient", col.levels = 8,
#' min.range = 0, max.range = 360)
#'
#' @export
dkcrop <- function(mesh, y){
  # Perform various checks:
  if (!isa(mesh, what = "mesh3d")) stop("mesh must be an object of class 'mesh3d'")
  if (!is.vector(y) & !is.integer(y) & !isa(y, what = "polygon.network")) stop("y should be a vector of integers or an object of class 'polygon.network'")
  if (is.vector(y)) {
    if(max(y) > Rvcg::nfaces(mesh)) stop(paste("maximum y (", max(y),") is larger than mesh face count (", Rvcg::nfaces(mesh), ")" ))
  }
  # Main job
  if (isa(y, what = "polygon.network")) y <- y@faces
  Croppedit <- sort(unique(c(mesh$it[, y])))
  Newvb <- mesh$vb[, Croppedit]
  dd <- data.frame(oldits = Croppedit, newits = c(1:length(Croppedit)))
  Newit <- mesh$it[, y]
  for (i in c(1:length(dd[,1]))){
    Newit[1, which(Newit[1,] == dd$oldits[i])] <- dd$newits[i]
    Newit[2, which(Newit[2,] == dd$oldits[i])] <- dd$newits[i]
    Newit[3, which(Newit[3,] == dd$oldits[i])] <- dd$newits[i]
  }
  NewMesh <- rgl::tmesh3d(vertices = Newvb, indices = Newit)
  NewMesh <- Rvcg::vcgUpdateNormals(NewMesh)
  #rgl::tmesh3d generates an object of class c('mesh3d', 'shape3d')
  #class('shape3d') needsto be removed for compatibility issues:
  class(NewMesh) <- "mesh3d"
  return(NewMesh)
}

# dkmap----
#' @title 3d topographic map
#' @description Map topographic variables on a 3d triangle mesh.
#' @param mesh object of class mesh3d
#' @param y a vector of values for each polygon of the mesh, usually a topographic variable
#' @param alpha an integer between 0 and 1 corresponding to alpha value for the selected polygons
#' (see \code{alpha.above}, \code{alpha.faces} and \code{alpha.thresh})
#' @param alpha.above logical, if TRUE polygons affected by \code{alpha} should have a y value above \code{alpha.thresh},
#' if FALSE their y value should be below \code{alpha.thresh}
#' @param alpha.faces a numeric vector of indices indicating which faces affected by \code{alpha}
#' @param alpha.thresh a numeric value indicating a threshold for alpha
#' @param bbox a logical, if TRUE a bounding box will be displayed around the surface object
#' @param bg the color to be used for the background. Defaults to "white".
#' @param cex a numerical value giving the amount by which plotting text and symbols should be magnified
#' relative to the default. This starts as 1 when a device is opened, and is reset when the layout is
#' changed, e.g. by setting mfrow.
#' @param cex.axis the magnification to be used for legend axis annotation relative to the current setting of cex.
#' @param cex.lab the magnification to be used for legend labels relative to the current setting of cex.
#' @param cex.main the magnification to be used for main titles relative to the current setting of cex.
#' @param cex.sub the magnification to be used for sub-titles relative to the current setting of cex.
#' @param col a vector of colors for texturing the polygons according to y
#' @param col.levels the number of color levels
#' @param col.axis the color to be used for legend axis annotation. Defaults to "black".
#' @param col.lab the color to be used for the legend labels. Defaults to "black".
#' @param col.main the color to be used for legend main titles. Defaults to "black".
#' @param col.sub the color to be used for plot sub-titles. Defaults to "black".
#' @param family the name of a font family for drawing text. The maximum allowed length is 200
#' bytes.
#' This name gets mapped by each graphics device to a device-specific font description.
#' The default value is "" which means that the default device fonts will be used (and what those
#' are should be listed on the help page for the device).
#' Standard values are "serif", "sans" and "mono", and the Hershey font families are also
#' available.
#' @param font.axis the font to be used for axis annotation.
#' @param font.lab the font to be used for the legend axis
#' @param font.main the font to be used for plot main titles.
#' @param font.sub the font to be used for plot sub-titles.
#' @param legend a logical indicating whether a legend should be displayed.
#' @param legend.lab a label for the legend axis.
#' @param legend.type a character string specifying the type of legend to be used; default is
#' "stack", which corresponds to a stacked vertical legend; "pie" generates a pie-shaped legend and "log" generates
#' a stacked vertical legend, but does a log transformation of the data (base: e=exp(1)). The "log" is mostly useful
#' for DNE maps.
#' @param lit logical, specifying if lighting calculation should take place on geometry
#' @param main the main title (on top) using font, size (character expansion) and color par(c("font.main", "cex.main", "col.main"))
#' @param sub sub-title (at bottom) using font, size and color par(c("font.sub", "cex.sub", "col.sub"))
#' @param max.range optional; the maximal range of the scale
#' @param min.range optional; the minimal range of the scale
#' @param orient the orientation of the view. For more details, see \code{\link{dksetview}}
#' @param origin logical, whether to set the z of the mesh's lowermost point to zero
#' @param scalebar A logical indicating whether a scalebar should be displayed
#' @param smooth A logical indicating whether the color of polygons should blend with neighbor polygons for a smoother rendering
#' @param windowRect the dimensions of the rgl window (default is the current size or, if size is below 1000*800, c(20, 20, 1020, 820))
#' @return An rgl window displaying the topography of a variable over a 3d mesh.
#' @seealso \code{\link[rgl]{rgl}}
#' @examples
#' #Map of orientation:
#' orient <- orient(dkmodel$complex)
#' dkmap(dkmodel$complex, orient, col.levels = 8, col = "orient",
#' legend.lab = "Orientation (degrees)",legend.type = "pie", min.range = 0,
#' max.range = 360, orient = "occlusal")
#'
#' #Map of area-relative curvature:
#' arc <- arc(dkmodel$complex)
#' dkmap(dkmodel$complex, arc, col = "arc", legend.lab = "ARC",
#' min.range = -20, max.range = 20, col.levels = 15, orient = "occlusal")
#'
#' #Map of Dirichlet normal energy:
#' dne <- dne(dkmodel$complex)
#' dkmap(dkmodel$complex, dne, col = "dne", legend.lab = "DNE",
#' legend.type = "log", orient = "occlusal")
#'
#' #Map of 3d-Area of polygons (for surface checking):
#' dkmap(dkmodel$complex, Rvcg::vcgArea(dkmodel$complex, perface = TRUE)$pertriangle,
#' legend.lab = "3d Area (mm\U00B2)", orient = "occlusal")
#' @export
dkmap <- function(mesh, y,  alpha = 1, alpha.above = TRUE, alpha.faces = NULL, alpha.thresh = NULL,
                   bg = "white", col = "slope", col.levels = 100, col.main = "black", col.lab = "black", col.sub = "black", col.axis = "black",
                   max.range = NULL, min.range = NULL, lit = TRUE,
                   cex = 2, cex.axis = 2, cex.main = 4, cex.sub = 3, cex.lab = 2,
                   family = "sans", font.axis = 1, font.lab = 2, font.main = 3, font.sub = 2, main = "", sub = "",
                   legend = TRUE, legend.lab = "y", legend.type = "stack", windowRect = NULL,
                   orient = "current", bbox = FALSE, origin = TRUE, scalebar = FALSE, smooth = TRUE){

  # Perform various checks:
  if (!isa(mesh, what = "mesh3d")) stop("mesh must be an object of class 'mesh3d'")
  if (!isa(y, what = "numeric") | !is.vector(y)) stop("y must be a numeric vector")
  if (!is.null(alpha.faces) & !is.null(alpha.thresh)) stop("use either alpha.faces or alpha.thresh")
  if (!is.null(alpha.faces) & (!isa(alpha.faces, what = "integer") | !is.vector(alpha.faces))) stop("alpha.faces must be a vector of integers")
  # Define colors
  Colrange <- col
  if (isTRUE(Colrange == "angularity")) Colrange <- c("white", "black")
  if (isTRUE(Colrange == "arc")) Colrange <- c("royalblue", "white", "red")
  if (isTRUE(Colrange == "dne")) Colrange <- c("royalblue", "lightskyblue", rep("olivedrab3", 3), "yellow1", "orange", "red")
  if (isTRUE(Colrange == "elev")) Colrange <- c("lightgreen","goldenrod1","yellow1","white","white","lightskyblue","dodgerblue4","royalblue")
  if (isTRUE(Colrange == "inclin")) Colrange <- c("firebrick4","red","orangered","orange","yellow1","olivedrab3","lightseagreen","royalblue","royalblue4","royalblue","lightseagreen","olivedrab3","yellow1","orange","orangered","red","firebrick4")
  if (isTRUE(Colrange == "oedist")) Colrange <- c("blue","green","yellow","orange","red")
  if (isTRUE(Colrange == "opc")) Colrange <- c("dodgerblue4","lightskyblue","sienna4","yellow1","red3","plum1","darkgreen","olivedrab3")
  if (isTRUE(Colrange == "opcr")) Colrange <- c("dodgerblue4","lightskyblue","sienna4","yellow1","red3","plum1","darkgreen","olivedrab3")
  if (isTRUE(Colrange == "orient")) Colrange <- c("dodgerblue4","lightskyblue","sienna4","yellow1","red3","plum1","darkgreen","olivedrab3")
  if (isTRUE(Colrange == "slope")) Colrange <- c("royalblue4","royalblue","lightseagreen","olivedrab3","yellow1","orange","orangered","red","firebrick4")

  # ...log transformation if legend type "log"
  if(legend.type == "log"){
    y <- log(y)
    #round -Inf and NA values to the lowest finite value
    y[!is.finite(y)] <- min(y[is.finite(y)])
  }

  # ...Define color palette of length = col.levels within col
  Colpalette <- grDevices::colorRampPalette(Colrange)(col.levels)
  # ...Define n = col.levels incremential levels for y
  y.min <- min(y)
  y.max <- max(y)
  if (!is.null(min.range)) y.min <- min.range
  if (!is.null(max.range)) y.max <- max.range
  Levels <- seq(y.min, y.max, length.out = col.levels + 1)
  # ...Assign a color of Colpalette to each polygon based on y level
  col.tmp <- c(1:length(y))
  for (i in c(1:col.levels)){
    col.tmp[y >= Levels[i]] <- Colpalette[i]
  }

  # Define alpha-ed faces
  # ...by default, to.alpha = no faces
  to.alpha <- c(1:Rvcg::nfaces(mesh))
  if (!is.null(alpha.faces)) to.alpha <- alpha.faces
  if (!is.null(alpha.thresh) & alpha.above == TRUE) to.alpha <- which(y > alpha.thresh)
  if (!is.null(alpha.thresh) & alpha.above == FALSE) to.alpha <- which(y < alpha.thresh)
  alpha.tmp <- rep(1, length(y))
  alpha.tmp[to.alpha] <- alpha

  # Populate the rgl window with the map
  # ...Clear the rgl window, resize it
  rgl::clear3d()
  if (!is.null(windowRect)) rgl::par3d(windowRect = windowRect)
  else if (rgl::par3d()$viewport[3] < 800 | rgl::par3d()$viewport[4] < 600) rgl::par3d(windowRect = c(20, 20, 1020, 820))
  # ...Set back origin to 0
  if (origin) {
    mesh <- doolkit::dkorigin(mesh)
  }
  # ...3dmap
  rgl::shade3d(mesh, alpha = alpha.tmp, color = col.tmp, meshColor = "faces", ptsize = 3, shininess = 0, smooth = smooth, lit = lit)

  if (lit) {
    # Light: 3 points lightning (https://fr.wikipedia.org/wiki/%C3%89clairage_trois_points)
    rgl::clear3d("lights") #sets the lights off
    ##key light
    rgl::light3d(specular = "black", diffuse = "grey50", ambient = "grey50", theta = 330, phi = 330)
    ##fill light
    rgl::light3d(specular = "black", diffuse = "grey60", ambient = "grey60", theta = 60, phi = 330)
    ##back light
    rgl::light3d(specular = "black", diffuse = "grey75", ambient = "grey50", theta = 315, phi = 45)
  }


  # Options
  # ...Legend
  if (legend == TRUE) {
    if (legend.type == "stack") {
      plotx <- c(3, rep(0, col.levels*100))
      ploty <- c(min(Levels), seq(min(Levels), max(Levels), length.out = col.levels*100))
      # ...function to make the bg
      bg.make <- function(){
        oldpar <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(oldpar))
        graphics::par(mar = c(3 + cex.sub, 3 + cex.lab + cex.axis, 3 + cex.main, 0), bg = bg)
        graphics::plot(x = plotx, y = ploty, pch = 95, bty = "n", xaxs = "r", xaxt = "n",
                       col = c("transparent", rep(Colpalette, each = 100)),
                       xlab = "", ylab = legend.lab, main = main, sub = sub,
                       col.main = col.main, col.axis = col.axis, col.lab = col.lab, col.sub = col.sub,
                       cex = cex, cex.main = cex.main, cex.axis = cex.axis, cex.sub = cex.sub, cex.lab = cex.lab,
                       family = family, font.main = font.main, font.axis = font.axis)
      }
      rgl::bgplot3d(bg.make())
    }

    if(legend.type == "log"){
      plotx <- c(3, rep(0, col.levels*100))
      ploty <- c(min(Levels), seq(min(Levels), max(Levels), length.out = col.levels*100))
      # ...function to make the bg
      bg.make <- function(){
        oldpar <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(oldpar))
        graphics::par(mar = c(3 + cex.sub, 3 + cex.lab + cex.axis, 3 + cex.main, 0), ylog = TRUE, bg = bg)
        graphics::plot(x = plotx, y = ploty, pch = 95, bty = "n", xaxs = "r", xaxt = "n",
                       col = c("transparent", rep(Colpalette, each = 100)),
                       xlab = "", ylab = paste("log(", legend.lab, ")", sep = ""), main = main, sub = sub,
                       col.main = col.main, col.axis = col.axis, col.lab = col.lab, col.sub = col.sub,
                       cex = cex, cex.main = cex.main, cex.axis = cex.axis, cex.sub = cex.sub, cex.lab = cex.lab,
                       family = family, font.main = font.main, font.axis = font.axis)
      }
      rgl::bgplot3d(bg.make())
    }

    if(legend.type == "pie"){
      ploty <- rep(1, col.levels)
      # ...function to make the bg
      bg.make <- function(){
        #... title
        oldpar <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(oldpar))
        graphics::par(mar = c(3 + cex.sub, 3 + (3 * cex), 3 + cex.main, 0), plt = c(0.15, 0.85, 0.15, 0.85), bg = bg)
        graphics::plot(1, 1, col = "white", main = main, font.main = font.main, cex.main = cex.main, col.main = col.main,
                       sub = sub, font.sub = font.sub, cex.sub = cex.sub, col.sub = col.sub,
                       family = family, xaxt = "n", yaxt = "n", bty = "n")
        #... piechart (just colors, without labels)
        graphics::par(plt = c(0.05, 0.2, 0.4, 0.6), new = TRUE)
        graphics::pie(ploty, labels = "", col = Colpalette, border = col.axis,
                      clockwise = TRUE, radius = 0.6)
        #... labels
        graphics::par(plt = c(0.05, 0.2, 0.4, 0.6), new = TRUE)
        graphics::pie(ploty, labels = round(seq(min(y), max(y), length.out = col.levels + 1), 0), init.angle = 90 + (360 / (2 * col.levels)),
                      clockwise = TRUE, radius = 0.8, col = "transparent", border = "transparent",
                      main = legend.lab, family = family, cex = cex, tck = NA,
                      font.main = font.lab, font.axis = font.axis,
                      cex.main = cex.lab, cex.axis = cex.axis,
                      col.main = col.lab, col.axis = col.axis)
      }
      rgl::bgplot3d(bg.make())
    }
  }

  else if (legend == FALSE) rgl::bgplot3d(NULL)

  #....3d view orientation
  if(orient == "occlusal" | orient == "mesial" | orient == "distal" | orient == "left" | orient == "right") {
    doolkit::dksetview(orient = orient)
  }
  if(orient != "current" & orient != "occlusal" & orient != "mesial" & orient != "distal" & orient != "left" & orient != "right") {
    print(paste("Warning: 'orient' argument must be of type c('current', 'distal', 'left', 'mesial', 'occlusal', 'right'); orient set to 'current'"))
  }

  #....Scalebar
  if(scalebar == TRUE) rgl::axes3d(c("x--", "y--", "z--"))
  else if(scalebar == FALSE)  rgl::axes3d(NULL)

  #....Bounding box
  if(bbox == TRUE) rgl::axes3d(edges = "bbox", labels = FALSE, tick = FALSE, box = TRUE, expand = 1, col = c("tan3", "white"))

  # The final step is setting the parallax to 0:
  rgl::par3d(FOV = 0)
}

# dkprofile----
#' @title cumulative profile, its slope and the area under its curve
#' @description A function for drawing the cumulative profile of a variable, computing the area under the curve and
#' the slope of the profile at the arithmetic mean of the variable.
#' @param x a numeric vector
#' @param type a character string indicating the type of coordinates to use
#' ("cartesian", "polar" etc.). Currently only "cartesian" is supported.
#' @param xlab title of the x axis
#' @param ylab title of the y axis
#' @param main main title of the plot
#' @param col the color of data points
#' @param alpha numeric indicating the alpha value of data points
#' @param size the size of data points
#' @param linetype the type of line to be traced (see ggplot2)
#' @return A list containing (1) the area under the curve of the profile, (2) the profile to be drawn,
#' and (3) the slope of the profile at the mean of the variable.
#' @references \doi{10.3389/fphys.2017.00524}{Thiery et al. (2017)}
#' @examples
#' #Elevation (hypsometric) profile (see Thiery et al., 2017):
#' dkprofile(elev(dkpongo$OES), main = "Elevation profile - Pongo pygmaeus",
#' ylab = "Elevation (%)", col = "#0072B2", linetype = "solid")
#'
#' #Enamel-dentine distance (pachymetric) profile:
#' dkprofile(oedist(dkpongo$OES, dkpongo$EDJ),
#' main = "Enamel-dentine distance profile - Pongo pygmaeus", ylab = "Distance (%)",
#' col = "#F0E442", linetype = "dashed")
#'
#' #Curvature (kurtometric) profile:
#' dkprofile(Rvcg::vcgCurve(dkpongo$OES)$meanitmax,
#' main = "Curvature profile - Pongo pygmaeus", ylab = "Curvature (%)",
#' col = "#D55E00", linetype = "dotted")
#'
#' @export
dkprofile <- function (x, type = 'cartesian', xlab = paste("cumulated frequency (%)"), ylab = '', main = '', col = "red", alpha = 1, size = 1, linetype = "solid") {
  if (min(x) < 0) x <- x + (abs(min(x))) #for kurtometric profiles
  if (min(x) > 0) x <- x - min(x) #for profiles of variables that do not start from 0
  Sorted <- x/(max(x)-min(x))*100
  Sorted <- Sorted[order(Sorted)]
  Cumulated <- c(1:length(Sorted))*100/length(Sorted)
  Df <- data.frame(Cumulated, Sorted)
  # ggplot of the profile:
  PPlot <- ggplot2::ggplot(Df, ggplot2::aes(x = Cumulated, y = Sorted)) +
    ggplot2::geom_line(colour = col, size = size, alpha = alpha, linetype = linetype) +
    ggplot2::ggtitle(main) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)
  if (type == 'cartesian') {
    PPlot <- PPlot + ggplot2::coord_cartesian()
  }
  else if (type == 'polar') {
    print(paste("Not implemented yet; type set to 'cartesian'"))
    PPlot <- PPlot + ggplot2::coord_cartesian()
  }
  else {
    print(paste("Warning: 'type' argument must be of type c('default','polar'); type set to 'default'"))
    PPlot <- PPlot + ggplot2::coord_cartesian()
  }
  # profile slope:
  Antecedent <- which.min(abs(Sorted - mean(Sorted)))
  Slopes <- NULL
  for (i in 1:round(length(x)/100, 0)) {
    Slopes[i] <- (Sorted[Antecedent + i] - Sorted[Antecedent - i])/(Cumulated[Antecedent + i] - Cumulated[Antecedent - i])
  }
  Slope <- round(mean(Slopes), 3)
  # area under curve:
  AUC <- sum(tis::lintegrate(Cumulated, Sorted, xint = Cumulated))

  return(list(auc = AUC, profile = PPlot, slope = Slope))
}

# dkorigin----
#' @title dkorigin
#' @description Sets the lowermost point of the mesh to 0 on the Z-axis
#' @param mesh object of class mesh3d
#' @return An object of class mesh3d.
#' @examples
#' #Map of elevation before using dkorigin:
#' dkmap(dkpongo$OES, doolkit::elev(dkpongo$OES), col = "elev", legend.lab = "Elevation (mm)")
#'
#' #Map of elevation after dkorigin:
#' leveled <- dkorigin(dkpongo$OES)
#' dkmap(leveled, doolkit::elev(leveled), col = "elev", legend.lab = "Elevation (mm)")
#' @export
dkorigin <- function(mesh){
  if (!isa(mesh, what = "mesh3d")) stop("mesh must be an object of class 'mesh3d'")
  Mesh <- mesh
  Mesh$vb[3, ] <- Mesh$vb[3, ] - min(mesh$vb[3, ])
  return(Mesh)
}

# dksetview----
#' @title preset orientations
#' @description A function to orient 3d topographical maps using preset values.
#' @param orient a character string indicating the targeted orientation (default is occlusal)
#' @return sets the orientation of the 'rgl' window.
#' @seealso \code{\link{dkmap}}
#' @examples
#' inclinCusp <- inclin(dkmodel$cusp)
#' dkmap(dkmodel$cusp, inclinCusp, col = "inclin", min.range = 0, max.range = 180)
#' dksetview()
#' #possible orientations are "distal", "left", "occlusal", "mesial" and "right"
#' @export
dksetview <- function(orient = "occlusal"){
  if (orient == "occlusal") rgl::view3d(theta = 0, phi = 0)
  if (orient == "mesial")   rgl::view3d(theta = 180, phi = 90)
  if (orient == "left")     rgl::view3d(theta = -90, phi = 0)
  if (orient == "right")    rgl::view3d(theta = 90, phi = 0)
  if (orient == "distal")   rgl::view3d(theta = 0, phi = -90)
  rgl::par3d(FOV = 0)
}
