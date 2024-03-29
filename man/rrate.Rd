% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topography.R
\name{rrate}
\alias{rrate}
\title{relief rate}
\usage{
rrate(uncropped, cropped, origin = TRUE)
}
\arguments{
\item{uncropped}{object of class mesh3d. Should entirely contain the 'cropped' argument.}

\item{cropped}{object of class mesh3d. Should be part of the 'uncropped' argument.}

\item{origin}{logical, if TRUE both cropped and uncropped mesh are translated along the z-axis
so that the lowest z of the uncropped mesh = 0; see \code{\link{dkorigin}}}
}
\value{
A single relief rate value.
}
\description{
Compute the relief rate from a sub-sample of a 3d triangle mesh.
For instance, the relief rate could be computed from the portion of a molar above
the lowermost point of its central basin, compared to the whole tooth.
}
\examples{
medelev <- median(elev(dkmodel$cusp))
basins <- dkcrop(dkmodel$cusp, which(elev(dkmodel$cusp) < medelev))
cusps <- dkcrop(dkmodel$cusp, which(elev(dkmodel$cusp) > medelev))

rrate(dkmodel$cusp, basins)
rrate(dkmodel$cusp, cusps)
}
\references{
\doi{10.1002/ajpa.23916}{Thiery et al. (2019)}
}
