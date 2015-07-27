\name{gDist}
\alias{gDist}
\alias{gDist.rtop}
\alias{gDist.SpatialPolygonsDataFrame}
\alias{gDist.SpatialPolygons}
\alias{gDist.list}
\title{calculate geostatistical distances between areas}
\description{Calculate geostatistical distances (Ghosh-distances) between
areas
}
\usage{
\method{gDist}{rtop}(object, params = list(), ...) 
\method{gDist}{SpatialPolygonsDataFrame}(object, object2 = NULL, ...) 
\method{gDist}{SpatialPolygons}(object, object2 = NULL, ...) 
\method{gDist}{list}(object, object2 = NULL, diag = FALSE, debug.level = 0, ...) 
}

\arguments{
\item{object}{object of class \code{\link[sp]{SpatialPolygons}} or 
\code{\link[sp:SpatialPolygons]{SpatialPolygonsDataFrame}} with boundaries of areas;
or list of discretized areas, typically from a call to \cr
\code{\link{rtopDisc}}; or
object of class \code{rtop} with such boundaries and/or discretized
elements (the individual areas)}
\item{params}{a set of parameters, used to modify the standard parameters for 
      the \code{rtop} package, set in \code{\link{getRtopParams}}. The argument params can also be used
      for the other methods, through the ...-argument.}
\item{object2}{an object of same type as \code{object}, except for \code{rtop};
for calculation of geostatistical
distances also between the elements in the two different objects }
\item{diag}{logical; if TRUE only calculate the geostatistical distances between
      each element and itself, only when the objects are lists of discretized
      areas and object2 = object or object2 = NULL}
\item{debug.level}{debug.level = 0 will suppress output from the call to 
      varMat, done for calculation of the geostatistical distances}
\item{...}{other parameters, for \code{gDist.list} when calling one of the 
other methods, or for \code{\link{varMat}}, in which the calculations take place}
} 

\value{ 
If called with one list of discretized elements, a matrix with the geostatistical distances 
between the elements within the list. If called with two lists of discretized elements, 
a matrix with the geostatistical distances between the elements in the two lists.
If called with \code{diag = TRUE}, the function returns an array of the geostatistical 
distance within each of the elements in the list.

If called with one \code{\link[sp]{SpatialPolygons}} or \code{\link[sp]{SpatialPolygonsDataFrame}}
or the function returns a list with one matrix with geostatistical distances between 
the elements of the object. If called with two objects, the list will also containt 
a matrix of the geostatistical distances between the elements of the two objects, and an array
of the geostatistical distances within the elements of the second object. 

If called with an rtop-object, the function will return the object, amended with the 
list above.
}

\note{
The geostatistical distance can be seen as the average distance between 
points in two elements, or the average distance within points in a single element.
The distance measure is also sometimes referred to as Ghosh-distance, from
Ghosh (1951) who found analytical expressions for these distances between 
blocks with regular geometry.

The use of geostatistical distances within \code{rtop} is based on an idea 
from Gottschalk (1993), who suggested
to replace the traditional regularization of variograms within block-kriging
(as done in the original top-kriging application of Skoien et al (2006))
with covariances of the geostatistical distance. The covariance between two
areas can then be found as \code{C(a1,a2) = cov(gd)} where \code{gd} is the geostatistical
distance between the two areas \code{a1} and \code{a2}, instead of an integration 
of the covariance function between the two areas. 

\code{rtop} is based on semivariograms
instead of covariances, and the semivariogram value between the two areas
can be found as \code{gamma(a1,a2) = g(gd) - 0.5 (g(gd1) + g(gd2))} where
\code{g} is a semivariogram valid for point support, \code{gd1)} and \code{gd2}
are the geostatistical distances within each of the two areas.

}


\references{ 

Ghosh, B. 1951. Random distances within a rectangle and between two rectangles. Bull. Calcutta Math. Soc., 43, 17-24.

Gottschalk, L. 1993. Correlation and covariance of runoff. Stochastic Hydrology and Hydraulics, 7, 85-101.

Skoien, J. O., R. Merz, and G. Bloschl. 2006. Top-kriging - geostatistics on stream networks. Hydrology and Earth System Sciences, 10, 277-287.


}
\author{ Jon Olav Skoien }
\examples{
\dontrun{
library(rgdal)
rpath = system.file("extdata",package="rtop")
observations = readOGR(rpath,"observations")
gDist = gDist(observations)
}

}
\keyword{spatial}
