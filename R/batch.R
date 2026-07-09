# batch.single
#' @title batch.single
#' @description Computes a range of dental topography variables for all triangles of a single mesh.
#' @param mesh An object of class mesh3d
#' @param functions A list of functions to apply; these functions should accept a 'mesh' argument and return as many values as the 'mesh' face count
#' @return A data frame.
#' @examples
#' ## Prepare mesh
#' my_mesh <- doolkit::dkmodel$cusp
#' ## Prepare functions
#' Slope <- function(mesh) return(doolkit::slope(mesh))
#' Rotation <- function(mesh) return(doolkit::orient(mesh))
#' fun_list <- list("slope" = Slope, "rotation" = Rotation)
#' ## Single mesh batch analysis
#' single_mesh_topography <- doolkit::batch.single(my_mesh, fun_list)
#' @export
batch.single <- function(mesh, functions){
  # Perform various checks:
  if (!isa(mesh, what = "mesh3d")) stop("mesh must be an object of class 'mesh3d'")
  for (fun in functions){
    if (!is.function(fun)) stop ("fun must be a valid method")
  }
  # Prepare dataset
  NFaces <- Rvcg::nfaces(mesh)
  Result <- data.frame("index" = c(1:NFaces))
  # Main loop
  for (fun in functions){
    FunResult <- fun(mesh)
    if (length(FunResult) < NFaces) stop ("fun must return as many values as the 'mesh' face count")
    Result <- cbind(Result, FunResult)
  }
  # Rename columns and return data frame
  colnames(Result) <- c("index", names(functions))
  return(Result)
  doolkit::dkpongo$OES
  Morpho::meshDist(doolkit::dkpongo$OES, doolkit::dkpongo$EDJ)
}

# batch.paired
#' @title batch.paired
#' @description Computes a range of dental topography variables for paired meshes.
#' @param meshA An object of class mesh3d; The mesh to which functions (ex. distance) are calculated.
#' @param meshB An object of class mesh3d; The mesh from which functions (ex. distance) are calculated.
#' @param functions A list of functions to apply; these functions should accept a 'mesh' argument and return as many values as the 'mesh' face count
#' @param method String indicating which method should be used to find triangle pairs.
#' @param direction Signed integer, indicates the direction for normals / orthogonal methods.
#' Default is '-1' which indicates that mesh B is below mesh A.
#' @param  epsilon Float corresponding to the error margin for intersections.
#' @return A data frame.
#' @examples
#' ## Prepare meshes
#' mesh_A <- doolkit::dkpongo$OES
#' mesh_B <- doolkit::dkpongo$EDJ
#' ## Prepare functions
#' Distance <- function(mesh) return(Morpho::meshDist(mesh_A, mesh_B, plot = FALSE))
#' SlopeDelta <- function(mesh) return(doolkit::slope(mesh_A) - doolkit::slope(mesh_B))
#' fun_list <- list("distance" = Distance, "slope_delta" = SlopeDelta)
#' ## Paired mesh batch analysis
#' paired_mesh_topography <- doolkit::batch.paired(mesh_A, mesh_B, "nearest", fun_list)
#' @export
batch.paired <- function(meshA, meshB, functions, method){
  # Perform various checks:
  if (!isa(meshA, what = "mesh3d")) stop("meshA must be an object of class 'mesh3d'")
  if (!isa(meshB, what = "mesh3d")) stop("meshB must be an object of class 'mesh3d'")
  for (fun in functions){
    if (!is.function(fun)) stop ("fun must be a valid method")
  }
  # Get tridX
  Result <- doolkit::tridx(meshA, meshB, method = method)
  NFaces <- length(Result[, 1])

  # Main loop
  for (fun in functions){
    FunResult <- fun(meshA, mesh_B_paired)
    if (length(FunResult) < NFaces) stop ("fun must return as many values as the 'meshA' face count")
    Result <- cbind(Result, FunResult)
  }
  # Rename columns and return data frame
  colnames(Result) <- c("mesh_A_index", "mesh_B_paired_index", names(functions))
  return(Result)
}


# batch.multi----
#' @title batch.multi
#' @description Computes a range of dental topography variables for a range of meshes.
#' @param files A list of surface files; accepted formats are as in RVCG::vcgImport
#' @param functions A list of functions to apply; these functions should accept a 'mesh' argument and return a single value
#' @return A data frame.
#' @seealso \code{\link[Rvcg]{vcgImport}}
#' @examples
#' ## Prepare a list of files
#' meshes_right <- list("pongo_OES" = doolkit::dkpongo$OES, "pongo_EDJ" = doolkit::dkpongo$EDJ)
#' ## Prepare functions
#' Mss <- function(mesh) return(mean(doolkit::slope(mesh)))
#' Rfi <- function(mesh) return(doolkit::rfi(mesh, method = "Boyer"))
#' Opcr <- function(mesh) return(doolkit::opcr(mesh)$opcr)
#' fun_list <- list("mss" = Mss, "rfi" = Rfi, "opcr" = Opcr)
#' ## Multiple meshes batch analysis
#' multi_mesh_topography <- doolkit::batch.multi(meshes_right, fun_list)
#' @export
#' @importFrom foreach %dopar%
batch.multi <- function(files, functions, filenames = NULL, do.parallel = TRUE){
  # Perform various checks:
  if (!is.null(filenames)){
    if (!is.vector(filenames)) stop("'filenames' must be a vector")
    else if (length(filenames) != length(files)) stop("'filenames' must have the same length as 'files'")
  }
  else {
    filenames <- sapply(strsplit(basename(files), "\\.(obj|ply|stl)$"), "[", 1)
  }

  for (i in 1:length(functions)){
    fun <- functions[[i]]
    if (!is.function(fun)) stop ("fun must be a valid method")
    TempResult <- fun(doolkit::dkmodel$cusp)
    if (is.null(names(functions)[i])) stop("all elements in the 'functions' list must be named")
    if (length(TempResult) != 1) stop ("functions passed to the 'functions' argument must return a single value")
  }

  # Prepare dataset
  Result <- data.frame()

  # Parallel
  if (do.parallel){
    cluster <- snow::makeSOCKcluster(parallel::detectCores() - 1)
    doSNOW::registerDoSNOW(cluster)
    nfiles <- length(files)
    file_count <- 0
    for (file in files){
      mesh <- Rvcg::vcgImport(file, silent = TRUE)
      if (!isa(mesh, what = "mesh3d")) stop("failed to convert an imported file to a 'mesh3d' object")
      Vector <- c(filenames[i], foreach::foreach(i = 1:length(functions), .combine = "c") %dopar% functions[[i]](mesh))
      Result <- rbind(Result, Vector)

      # Check if a shiny progress function has been registered in the background
      shiny_updater <- getOption("doolkit.progress_callback")
      if (is.function(shiny_updater)) {
        file_count <- file_count + 1
        shiny_updater(amount = 1/nfiles, text = paste("Step", file_count, "of", nfiles))
      }

    }
    parallel::stopCluster(cluster)
  }
  # Not parallel
  else {
    for (file in files){
      mesh <- Rvcg::vcgImport(file, silent = TRUE)
      if (!isa(mesh, what = "mesh3d")) stop("failed to convert an imported file to a 'mesh3d' object")
      Vector <- vector()
      for (fun in functions) {
        FunResult <- fun(mesh)
        Vector <- c(Vector, fun(mesh))
      }
      Result <- rbind(Result, Vector)
    }
  }

  # Rename
  colnames(Result) <- c("id", names(functions))
  if (!is.null(filenames)){
    Result[, 1] <- c(filenames)
  }

  return(Result)
}


# dksave----
#' @title dksave
#' @description Save a data frame.
#' @param df An object of class data.frame
#' @param file The path to which the data.frame must be saved
#' @param format The format to which the data frame should be saved, either ".txt", ".csv", ".RDS" or ".xlsx" (default is "txt"))
#' @return A 'format' file saved at the 'path' argument.
#' @examples
#' # Build a topographic data frame
#' mesh_right <- doolkit::dkpongo$OES
#' Slope <- function(mesh) return(doolkit::slope(mesh))
#' Rotation <- function(mesh) return(doolkit::orient(mesh))
#' fun_list <- list("slope" = Slope, "rotation" = Rotation)
#' tooth_topo_dataframe <- doolkit::tooth_topography(mesh_right, fun_list)
#' #
#' # Save the dataframe to file...
#' # ...as .txt:
#' tmpfile <- tempfile("topographic_table", fileext = ".txt")
#' dksave(tooth_topo_dataframe, tmpfile, format = ".txt")
#' # ...as .csv:
#' tmpfile <- tempfile("topographic_table", fileext = ".csv")
#' dksave(tooth_topo_dataframe, tmpfile, format = ".csv")
#' # ...as .RDS:
#' tmpfile <- tempfile("topographic_table", fileext = ".RDS")
#' dksave(tooth_topo_dataframe, tmpfile, format = ".RDS")
#' # ...as .xlsx:
#' tmpfile <- tempfile("topographic_table", fileext = ".xlsx")
#' dksave(tooth_topo_dataframe, tmpfile, format = ".xlsx")
#'
#' unlink(tmpfile)
#' @export
dksave <- function(df, file, format = ".txt"){
  # Perform various checks
  if (!isa(df, what = "data.frame")) stop("all objects in 'meshes' must be of class of class 'mesh3d'")
  if (!(format %in% c(".txt", ".csv", ".RDS", ".xlsx"))) stop(paste(format, " is not a valid format"))
  # save
  if (".txt" %in% format){
    write.table(df, file = paste0(file , ".txt"), sep = "\t", row.names = F, append = F, quote = F)
  }
  else if (".csv" %in% format){
    write.table(df, file = paste0(file , ".csv"), sep = "\t", row.names = F, append = F, quote = F)
  }
  else if (".RDS" %in% format){
    saveRDS(df, paste0(file , ".RDS"))
  }
  else if (".xlsx" %in% format){
    openxlsx::write.xlsx(df, file = paste0(file , ".xlsx"), asTable = T)
  }
}

