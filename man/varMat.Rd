\name{varMat}
\alias{varMat}
\alias{varMat.rtop}
\alias{varMat.SpatialPolygonsDataFrame}
\alias{varMat.SpatialPolygons}
\alias{varMat.matrix}
\alias{varMat.list}
\title{ create a semivariogram matrix between a set of locations, or semivariogram 
matrices between and within two sets of locations}
\description{
varMat will create a semivariogram matrix between all the supports in a set of 
locations (observations or prediction locations) or semivariogram matrices
between all the supports in one or two sets of locations, and also between them.
}
\usage{
\method{varMat}{rtop}(object, varMatUpdate = FALSE, fullPred = FALSE, params = list(), ...) 
\method{varMat}{SpatialPolygonsDataFrame}(object, object2 = NULL,...) 
\method{varMat}{SpatialPolygons}(object, object2 = NULL, variogramModel,
     overlapObs, overlapPredObs, ...) 
%\method{varMat}{matrix}(object, variogramModel,  ...) 
\method{varMat}{list}(object, object2 = NULL, coor1, coor2, maxdist = Inf, 
              variogramModel, diag = FALSE, sub1, sub2, 
              debug.level = ifelse(interactive(), 1, 0), ...) 
}


\arguments{
\item{object}{ either: 1) an object of class \code{rtop} (see \link{rtop-package}) 
         or 2) a \cr
         \code{\link[sp:SpatialPolygons]{SpatialPolygonsDataFrame}},
               or \code{\link[sp:SpatialPolygons]{SpatialPolygons}}, 
         or 3) a \cr
         \code{\link{matrix}} with geostatistical distances (see
               \code{\link{gDist}}
         or 4) a \code{\link{list}} with discretized supports}
\item{varMatUpdate}{logical; if TRUE, also existing variance matrices will 
            be recomputed, if FALSE, only missing variance matrices will be computed}
\item{fullPred}{logical; whether to create the full covariance matrix also for the predictions, 
                mainly used for simulations}
\item{params}{a set of parameters, used to modify the default parameters for 
      the \code{rtop} package, set in \code{\link{getRtopParams}}. }
\item{object2}{if \code{object} is not an object of class \code{rtop}; 
 an object of the same class as \code{object} with a possible second 
 set of locations with support}
\item{variogramModel}{variogramModel to be used in calculation of the 
   semivariogram matrix (matrices)}
\item{...}{typical parameters to modify from the default parameters of the
    rtop-package (or modifications of the previously set parameters for the 
    \code{rtop}-object), see also \code{\link{getRtopParams}}. These can also 
    be passed in a list named params, as for the rtop-method. Typical 
    parameters to modify for this function:
    \itemize{
      \item{rresol = 100}{miminum number of discretization points, in call to 
        \code{\link{rtopDisc}} if necessary}
      \item{rstype = "rtop"}{sampling type from areas, in call to 
        \code{\link{rtopDisc}} if necessary}
     \item{gDistPred = FALSE}{use geostatistical distance for semivariogram matrices}
     \item{gDist}{parameter to set jointly \code{gDistEst = gDistPred = gDist}}
    }
}
\item{overlapObs}{matrix with observations that overlap each other}
\item{overlapPredObs}{matrix with \code{observations} and \code{predictionLocations}
          that overlap each other}
\item{coor1}{coordinates of centroids of \code{object}}
\item{coor2}{coordinates of centre-of-gravity of \code{object2}}
\item{maxdist}{maximum distance between areas for inclusion in semivariogrma matrix}
\item{diag}{logical; if TRUE only the semivariogram values along the diagonal 
       will be calculated, typical for semivariogram matrix of prediction locations}
\item{sub1}{semivariogram array for subtraction of inner variances of areas}
\item{sub2}{semivariogram array for subtraction of inner variances of areas}
\item{debug.level}{debug.level >= 1 will give output for every element}
} 

\value{ 
The lower level versions of the function calculates a semivariogram matrix
between locations in \code{object} or between the locations in \code{object} 
and the locations in \code{object2}. The method for object of type \code{rtop}
calculates semivariogram matrices between observation locations, between prediction locations, 
and between observation locations and prediction locations, and adds these 
to \code{object}.

}


\note{
The argument \code{varMatUpdate} is typically used to avoid repeated computations
of the same variance matrices. Default is FALSE, which will avoid recomputation
of the variance matrix for the observations if the procedure is cross-validation
before interpolation. Should be set to TRUE if the variogram Model has been
changed, or if observation and/or prediction locations have been changed.

If an \code{rtop}-object contains observations and/or predictionLocations of type 
\code{\link[spacetime]{STSDF}}, the covariance matrix is computed based on the spatial 
properties of the object.
}

\references{ 
Skoien J. O., R. Merz, and G. Bloschl. Top-kriging - geostatistics on stream networks.
Hydrology and Earth System Sciences, 10:277-287, 2006.

Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione, A., 2014. Rtop: An R package for interpolation of data with a variable spatial support, with an example from river networks. Computers & Geosciences, 67.
}


\author{ Jon Olav Skoien }
\seealso{
\code{\link{gDist}}
}
\examples{
\dontrun{
library(rgdal)
rpath = system.file("extdata",package="rtop")
observations = readOGR(rpath,"observations")
vmod = list(model = "Ex1", params = c(0.00001,0.007,350000,0.9,1000))
vm = varMat(observations, variogramModel = vmod)
}
}
\keyword{spatial}
