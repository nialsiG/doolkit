% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{dkorigin}
\alias{dkorigin}
\title{dkorigin}
\usage{
dkorigin(mesh)
}
\arguments{
\item{mesh}{object of class mesh3d}
}
\value{
An object of class mesh3d.
}
\description{
Sets the lowermost point of the mesh to 0 on the Z-axis
}
\examples{
#Map of elevation before using dkorigin:
dkmap(dkpongo$OES, doolkit::elev(dkpongo$OES), col = "elev", legend.lab = "Elevation (mm)")

#Map of elevation after dkorigin:
leveled <- dkorigin(dkpongo$OES)
dkmap(leveled, doolkit::elev(leveled), col = "elev", legend.lab = "Elevation (mm)")
}
