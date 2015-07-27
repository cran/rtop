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
