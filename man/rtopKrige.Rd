\name{rtopKrige}
\alias{rtopKrige}
\alias{rtopKrige.rtop}
\alias{rtopKrige.SpatialPolygonsDataFrame}
\alias{rtopKrige.default}
\title{Spatial interpolation of data with spatial support}
\description{
rtopKrige perform spatial interpolation or cross validation of data with areal support.
}
\usage{
\method{rtopKrige}{rtop}(object, varMatUpdate = FALSE, ...) 
\method{rtopKrige}{SpatialPolygonsDataFrame}(object, predictionLocations = NULL,
    varMatObs, varMatPredObs, varMat, params = list(), formulaString, 
    sel, ...) 
\method{rtopKrige}{default}(object, predictionLocations = NULL,
    varMatObs,  varMatPredObs, varMat, params = list(), formulaString, 
    sel, wret = FALSE, ...) 
}

\arguments{
\item{object}{object of class \code{rtop} or \code{\link[sp]{SpatialPolygonsDataFrame}}}
\item{varMatUpdate}{logical; if TRUE, also existing variance matrices will 
            be recomputed, if FALSE, only missing variance matrices will be computed,
            see also \code{\link{varMat}}}
\item{predictionLocations}{\code{\link[sp]{SpatialPolygons}} or \code{\link[sp]{SpatialPolygonsDataFrame}}
     with prediction locations. NULL if cross validation is to be performed.}
\item{varMatObs}{covariance matrix of observations, where diagonal must consist 
      of internal variance, typically generated from call
      to \code{\link{varMat}} }
\item{varMatPredObs}{covariance matrix between observation locations and 
      prediction locations, typically generated from call
      to \code{\link{varMat}}}
\item{varMat}{list covariance matrices including the two above}
\item{params}{a set of parameters, used to modify the standard parameters for 
      the \code{rtop} package, set in \code{\link{getRtopParams}}. Additionally,
      it is possible overrule some of the parameters in \code{object$params} by passing
      them as separate arguments.}
\item{formulaString}{formula that defines the dependent variable as a linear model 
      of independent variables, see e.g. \code{\link{createRtopObject}} for more details.}
\item{sel}{array of prediction location numbers, if only a limited number of locations are to be
     interpolated/crossvalidated}
\item{wret}{logical; if TRUE, return a matrix of weights instead of the predictions,
useful for batch processing of time series, see also details}
\item{...}{from \code{rtopKrige.rtop}, arguments to be passed to 
    \code{rtopKrige.default}. In \code{rtopKrige.default},
    parameters for modification of the object parameters or default parameters.
    Of particular interest are \code{cv}, a logical for doing cross-validation,
    \code{nmax}, and \code{maxdist} for maximum number of neighbours and
    maximum distance to neighbours, respectively, and \code{wlim}, the limit for
    the absolute values of the weights.}
} 

\value{ 
If called with \code{\link[sp]{SpatialPolygonsDataFrame}}, the function returns 
a \code{\link[sp]{SpatialPolygonsDataFrame}} with predictions, either at the 
locations defined in \code{predictionLocations}, or as leave-one-out
cross-validation predicitons at the same locations as in object if 
\code{cv = TRUE}

If called with an rtop-object, the function returns the same object with the
predictions added to the object.
}

\details{
This function is the interpolation routine of the rtop-package.
The simplest way of calling the function is with an rtop-object that 
contains the fitted variogram model and all the other necessary data (see
\code{\link{createRtopObject}} or \code{\link{rtop-package}}). 

The function will, if called with covariance matrices between observations
and between observations and prediction locations, use these for the interpolation.
If the function is called without these matrices, \code{\link{varMat}} will be
called to create them. These matrices can therefore be reused if necessary,
an advantage as it is computationally expensive to create them.

The interpolation that takes part within \code{rtopKrige.default} is based on 
the semivariance matrices between observations and between observations and prediction
locations. It is therefore possible to use this function also to interpolate
data where the matrices have been created in other ways, e.g. based on distances
in physiographical space or distances along a stream.

The function returns the weights rather than the predictions if \code{wret = TRUE}.
This is useful for batch processing of time series, e.g. if 

}




\references{ 
Skoien J. O., R. Merz, and G. Bloschl. Top-kriging - geostatistics on stream networks.
Hydrology and Earth System Sciences, 10:277-287, 2006.
}
\author{ Jon Olav Skoien }
\seealso{\code{\link{rtop-package}}}
\examples{
\dontrun{
library(rgdal)
rpath = system.file("extdata",package="rtop")
observations = readOGR(rpath,"observations")
predictionLocations = readOGR(rpath,"predictionLocations")

# Setting some parameters 
params = list(gDist = TRUE, cloud = FALSE)
# Create a column with the specific runoff:
observations$obs = observations$QSUMMER_OB/observations$AREASQKM
# Build an object
rtopObj = createRtopObject(observations,predictionLocations, params = params)
# Fit a variogram (function also creates it)
rtopObj = rtopFitVariogram(rtopObj)
# Predicting at prediction locations
rtopObj = rtopKrige(rtopObj)
# Cross-validation
rtopObj = rtopKrige(rtopObj,cv=TRUE)
cor(rtopObj$predictions$observed,rtopObj$predictions$var1.pred)
}
}
\keyword{spatial}
