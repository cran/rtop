\name{createRtopObject}
\alias{createRtopObject}
\title{Create an object for interpolation within the rtop package}
\description{
This is a help function for creating an object (see \code{\link{rtop-package}}
 to be used for interpolation within the rtop package
}
\usage{
createRtopObject(observations, predictionLocations,
   formulaString, params=list(), overlapObs,
   overlapPredObs, ...)
}

\arguments{
\item{observations}{\code{\link[sp]{SpatialPolygonsDataFrame}} or \code{\link[sf]{sf}}-polygons with observations}
\item{predictionLocations}{a \code{\link[sp]{SpatialPolygons}},  
      \code{\link[sp]{SpatialPolygonsDataFrame}}-object or 
       \code{\link[sf]{sf}}-polygons with 
                 prediction locations}
\item{formulaString}{formula that defines the dependent variable as a linear model 
      of independent variables; suppose the dependent variable has name \code{z}, 
      for ordinary and simple kriging use the formula \code{z~1}; 
      for universal kriging, suppose \code{z} is linearly dependent on 
      \code{x} and \code{y}, use the formula \code{z~x+y}. The formulaString defaults
      to \code{"value~1"} if \code{value} is a part of the data set. 
      If not, the first column of the data set is used. Universal kriging 
      is not yet properly implemented in the \code{rtop}-package, this 
      element is mainly used for defining the dependent variable.}
\item{params}{parameters to modify the default parameters of the rtop-package, 
        set internally in this function by a call to \code{\link{getRtopParams}} }

\item{overlapObs}{matrix with observations that overlap each other}
\item{overlapPredObs}{matrix with \code{observations} and \code{predictionLocations}
          that overlap each other}
\item{...}{Extra parameters to \code{\link{getRtopParams}} and possibility
       to pass depreceted arguments}
}

\value{ 
An object of class \code{rtop} with observations, prediction locations, 
parameters and possible other elements useful for interpolation in the rtop-package.
Most other externally visible functions in the package will be able to 
work with this object, and add the results as a new element.
}


\references{ 
Skoien J. O., R. Merz, and G. Bloschl. Top-kriging - geostatistics on stream networks.
Hydrology and Earth System Sciences, 10:277-287, 2006.

Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione, A., 2014. Rtop: An R package for interpolation of data with a variable spatial support, with an example from river networks. Computers & Geosciences, 67.
}


\author{ Jon Olav Skoien }


\examples{
\donttest{
rpath = system.file("extdata",package="rtop")
# As rgdal is about to be retired:
if (require("rgdal")) {
  library(rgdal)
  observations = readOGR(rpath,"observations")
  predictionLocations = readOGR(rpath,"predictionLocations")
} else {
  library(sf)
  observations = st_read(rpath, "observations")
  predictionLocations = st_read(rpath,"predictionLocations")
}
# Create a column with the specific runoff:
observations$obs = observations$QSUMMER_OB/observations$AREASQKM

# Setting some parameters 
params = list(gDist = TRUE, cloud = FALSE)
# Create a column with the specific runoff:
observations$obs = observations$QSUMMER_OB/observations$AREASQKM
# Build an object
rtopObj = createRtopObject(observations, predictionLocations, 
                           params = params)
}

}

\seealso{
\code{\link{getRtopParams}}
}
\keyword{spatial}
