# dkpongo----
#' @title dkpongo
#' @description A dataset containing the OES and the EDJ surfaces of a \emph{Pongo pygmaeus} upper second molar (SMF-1117)
#' @example Pongo_OES <- dkpongo$OES
#' Pongo_EDJ <- dkpongo$EDJ
#' @source \url{https://www.morphosource.org/Detail/MediaDetail/Show/media_id/42357}
dkpongo <- list(OES = Rvcg::vcgPlyRead(system.file("extdata", "plys", "Pongo_pygmaeus_SMF-1117_EM_50k.ply", package = "doolkit"), updateNormals = TRUE, clean = TRUE),
                EDJ = Rvcg::vcgPlyRead(system.file("extdata", "plys", "Pongo_pygmaeus_SMF-1117_DE_50k.ply", package = "doolkit"), updateNormals = TRUE, clean = TRUE))
usethis::use_data(dkpongo, overwrite = TRUE)
'dkpongo'

# dkmodel----
#' @title dkmodel
#' @description A list containing theoretical model surfaces corresponding to a flat surface, a single cusp, a shallow basin and a complex surface.
#' @example Flat_surface <- dkmodel$flat
#' Single_cusp <- dkmodel$cusp
#' Shallow_basin <- dkmodel$basin
#' Complex_surface <- dkmodel$complex
#' @source \url{https://github.com/nialsiG/A-comparison-of-relief-estimates-used-in-3d-dental-topography}
dkmodel <- list(basin = Rvcg::vcgPlyRead(system.file("extdata", "plys", "basin1.ply", package = "doolkit"), updateNormals = TRUE, clean = TRUE),
                complex = Rvcg::vcgPlyRead(system.file("extdata", "plys", "complex.ply", package = "doolkit"), updateNormals = TRUE, clean = TRUE),
                cusp = Rvcg::vcgPlyRead(system.file("extdata", "plys", "cusp.ply", package = "doolkit"), updateNormals = TRUE, clean = TRUE),
                flat = Rvcg::vcgPlyRead(system.file("extdata", "plys", "flat.ply", package = "doolkit"), updateNormals = TRUE, clean = TRUE))
usethis::use_data(dkmodel, overwrite = TRUE)
'dkmodel'
