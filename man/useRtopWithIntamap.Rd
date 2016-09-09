\name{useRtopWithIntamap}
\alias{useRtopWithIntamap}
\title{
Integrates the rtop package with the \code{intamap} package
}
\description{
This function makes it possible to use \code{rtop}-objects in the functions of the package.
It is necessary to load the \code{intamap}-package before calling this function. 
}
\usage{
useRtopWithIntamap()
}

\references{ 
Pebesma, E., Cornford, D., Dubois, G., Heuvelink, G.B.M., Hristopulos, D., Pilz, J., Stohlker, U., Morin, G., Skoien, J.O. INTAMAP: The design and implementation f an interoperable automated interpolation Web Service. Computers and Geosciences 37 (3), 2011. 

Skoien J. O., R. Merz, and G. Bloschl. Top-kriging - geostatistics on stream networks.
Hydrology and Earth System Sciences, 10:277-287, 2006.

Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione, A., 2014. Rtop: An R package for interpolation of data with a variable spatial support, with an example from river networks. Computers & Geosciences, 67.
}


\value{
The function will have as side effect that the intamap package is loaded, 
and that rtop-methods are registered for the intamap-functions estimateParameters,
spatialPredict and methodParameters.
}


\author{ Jon Olav Skoien }
\examples{
library(intamap)
useRtopWithIntamap()
}
\keyword{plot}
