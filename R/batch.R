# tooth_topography
#' @title tooth_topography
#' @description Computes a range of dental topography variables for all triangles of a single mesh.
#' @param mesh An object of class mesh3d
#' @param functions A list of functions to apply; these functions should accept a 'mesh' argument and return as many values as the 'mesh' face count
#' @return A data frame.
#' @examples
#' ## Prepare mesh
#' mesh_right <- doolkit::dkpongo$OES
#' ## Prepare functions
#' Slope <- function(mesh) return(doolkit::slope(mesh))
#' Rotation <- function(mesh) return(doolkit::orient(mesh))
#' fun_list <- list("slope" = Slope, "rotation" = Rotation)
#' ## Batch tooth analysis
#' tooth_topo_dataframe <- doolkit::tooth_topography(mesh_right, fun_list)
#' @export
tooth_topography <- function(mesh, functions){
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
}

# sample_topography----
#' @title sample_topography
#' @description Computes a range of dental topography variables for a range of meshes.
#' @param meshes A list of objects of class mesh3d
#' @param functions A list of functions to apply; these functions should accept a 'mesh' argument and return a single value
#' @return A data frame.
#' @examples
#' ## Prepare meshes
#' meshes_right <- list("pongo_OES" = doolkit::dkpongo$OES, "pongo_EDJ" = doolkit::dkpongo$EDJ)
#' ## Prepare functions
#' Mss <- function(mesh) return(mean(doolkit::slope(mesh)))
#' Rfi <- function(mesh) return(doolkit::rfi(mesh, method = "Boyer"))
#' Opcr <- function(mesh) return(doolkit::opcr(mesh)$opcr)
#' fun_list <- list("mss" = Mss, "rfi" = Rfi, "opcr" = Opcr)
#' ## Batch sample analysis
#' sample_topo_dataframe <- doolkit::sample_topography(meshes_right, fun_list)

#' @export
sample_topography <- function(meshes, functions){
  # Perform various checks:
  for (mesh in meshes){
    if (!isa(mesh, what = "mesh3d")) stop("all objects in 'meshes' must be of class of class 'mesh3d'")
  }
  for (fun in functions){
    if (!is.function(fun)) stop ("fun must be a valid method")
  }
  # Prepare dataset
  Result = data.frame("mesh" = names(meshes))
  # Main loop
  for (fun in functions){
    Vector <- vector()
    for (mesh in meshes){
      FunResult <- fun(mesh)
      if (length(FunResult) != 1) stop ("fun must return a single value")
      Vector <- c(Vector, FunResult)
    }
    Result <- cbind(Result, Vector)
  }
  # Rename columns and return data frame
  colnames(Result) <- c("mesh", names(functions))
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

