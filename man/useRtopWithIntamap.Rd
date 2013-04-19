\name{useRtopWithIntamap}
\alias{useRtopWithIntamap}
\title{
Integrates the rtop package with the intamap package
}
\description{
This function checks if the \code{\link[intamap]{intamap-package}} is installed,
and if it is, makes it possible to use \code{rtop}-objects in the functions of the package.
}
\usage{
useRtopWithIntamap()
}

\value{
The function will have as side effect that the intamap package is loaded, 
and that rtop-methods are registered for the intamap-functions estimateParameters,
spatialPredict and methodParameters.
}


\author{ Jon Olav Skoien }
\examples{
useRtopWithIntamap()
}
\keyword{plot}
