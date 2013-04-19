\name{createRtopObject}
\alias{createRtopObject}
\title{Create an object for interpolation within the rtop package}
\description{
This is a help function for creating an object (see \code{\link{rtop-package}}
 to be used for interpolation within the rtop package
}
\usage{
createRtopObject(observations, predictionLocations,
   formulaString, params=list(), ainfo, areas, overlapObs, overlapPredObs, ...)
}

\arguments{
\item{observations}{\code{\link[sp]{SpatialPolygonsDataFrame}} with observations}
\item{predictionLocations}{a \code{\link[sp]{SpatialPolygons}} or a 
      \code{\link[sp]{SpatialPolygonsDataFrame}}-object with 
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
\item{params}{parameters to modify the standard parameters of the rtop-package, 
        set internally in this function by a call to \code{\link{getRtopParams}} }

\item{ainfo}{\code{\link[sp]{SpatialPointsDataFrame}} with information about the 
      observations. Only used here if it contains information for separating
      the observations into observations and predictionLocations or
      combining \code{ainfo} with \code{areas} to find \code{observations} and 
      \code{predictionLocations}(depreceted method). Typically from a call to
      \code{\link{readAreaInfo}} }
\item{areas}{\code{\link[sp]{SpatialPolygonsDataFrame}} with areal information,
      typically from a call to \code{\link{readAreas}}, depreceted}
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



\author{ Jon Olav Skoien }

\seealso{
\code{\link{rtop-package}} and \code{\link{getRtopParams}}
}
\keyword{spatial}
