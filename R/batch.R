# tooth_topography
#' @title tooth_topography
#' @description Computes a range of dental topography variables for all triangles of a single mesh.
#' @param mesh An object of class mesh3d
#' @param functions A list of functions to apply; these functions should accept a 'mesh' argument and return as many values as the 'mesh' face count
#' @return A data frame.
#' @examples
#' #
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

  # Rename columns
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
#' #
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

  # Rename columns
  colnames(Result) <- c("mesh", names(functions))

  return(Result)

}
