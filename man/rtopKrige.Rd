\name{rtopKrige}
\alias{rtopKrige}
\alias{rtopKrige.rtop}
\alias{rtopKrige.SpatialPolygonsDataFrame}
\alias{rtopKrige.STSDF}
\alias{rtopKrige.default}
\title{Spatial interpolation of data with spatial support}
\description{
rtopKrige perform spatial interpolation or cross validation of data with areal support.
}
\usage{
\method{rtopKrige}{rtop}(object, varMatUpdate = FALSE, params = list(), ...) 
\method{rtopKrige}{SpatialPolygonsDataFrame}(object, predictionLocations = NULL,
    varMatObs, varMatPredObs, varMat, params = list(), 
    formulaString, sel, ...) 
\method{rtopKrige}{STSDF}(object, predictionLocations = NULL,
    varMatObs, varMatPredObs, varMat, params = list(), 
    formulaString, sel, olags = NULL, plags = NULL, 
    lagExact = TRUE, ...) 
\method{rtopKrige}{default}(object, predictionLocations = NULL,
    varMatObs,  varMatPredObs, varMat, params = list(),
    formulaString, sel, wret = FALSE, ...) 



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
\item{params}{a set of parameters, used to modify the default parameters for 
      the \code{rtop} package, set in \code{\link{getRtopParams}}. Additionally,
      it is possible overrule some of the parameters in \code{object$params} by passing
      them as separate arguments.}
\item{formulaString}{formula that defines the dependent variable as a linear model 
      of independent variables, see e.g. \code{\link{createRtopObject}} for more details.}
\item{sel}{array of prediction location numbers, if only a limited number of locations are to be
     interpolated/crossvalidated}
\item{wret}{logical; if TRUE, return a matrix of weights instead of the predictions,
useful for batch processing of time series, see also details}
\item{olags}{A vector describing the relative lag which should be applied for the observation locations. See also details}
\item{plags}{A vector describing the relative lag which should be applied for the predicitonLocations. See also details}
\item{lagExact}{logical; whether differences in lagtime should be computed exactly or approximate}
\item{...}{from \code{rtopKrige.rtop}, arguments to be passed to 
    \code{rtopKrige.default}. In \code{rtopKrige.default},
    parameters for modification of the object parameters or default parameters.
    Of particular interest are \code{cv}, a logical for doing cross-validation,
    \code{nmax}, and \code{maxdist} for maximum number of neighbours and
    maximum distance to neighbours, respectively, and \code{wlim}, the limit for
    the absolute values of the weights.}
} 

\value{ 
If called with \code{\link[sp]{SpatialPolygonsDataFrame}}, the function returns a \cr
\code{\link[sp]{SpatialPolygonsDataFrame}} with predictions, either at the 
locations defined in \cr
\code{predictionLocations}, or as leave-one-out
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
This is useful for batch processing of time series, e.g. once the weights are 
created, they can be used to compute the interpolated values for each time step.

rtop is able to take some advantage of multiple cpus, which can be invoked with the 
parameter nclus, setting this larger than one will start a cluster with nclus workers,
if the \code{\link{parallel}}-package has been installed. 


Kriging of time series is possible when \code{observations} and \code{predictionLocations}
are spatiotemporal objects of type \code{\link[spacetime]{STSDF}}. The interpolation is 
still spatial, in the sense that the regularisation of the variograms are just done
using the spatial extent of the observations, not a possible temporal extent, such as
done by Skoien and Bloschl (2007). However, it is possible to make predictions based on observations 
from different time steps, through the use of the lag-vectors. These vectors describe a typical "delay"
for each observation and prediction location. This delay could for runoff related variables be similar
to travel time to each gauging location. For a certain prediction location, earlier time steps would be picked for neighbours with shorter travel time and later time steps for neighbours with slower travel times.

The lagExact parameter indicates whether to use a weighted average of two time steps, or just the time step which is closest to the difference in lag times. 

The use of lag times should in theory increase the computation time, but might, due to different computation methods, even speed up the computation when the number of neighbours to be used (parameter nmax) is small compared to the number of observations. If computation is slow, it can be useful to test olags = rep(0, dim(observations)[1]) and similar for predictionLocations.




}




\references{ 
Skoien J. O., R. Merz, and G. Bloschl. Top-kriging - geostatistics on stream networks.
Hydrology and Earth System Sciences, 10:277-287, 2006.

Skoien, J. O. and G. Bloschl. Spatio-Temporal Top-Kriging of Runoff Time Series. Water Resources Research 43:W09419, 2007.

Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione, A., 2014. Rtop: An R package for interpolation of data with a variable spatial support, with an example from river networks. Computers & Geosciences, 67.
}

\author{ Jon Olav Skoien }
\examples{
\dontrun{
library(rgdal)
# The following command will download  the complete example data set
# downloadRtopExampleData() 
# observations$obs = observations$QSUMMER_OB/observations$AREASQKM

rpath = system.file("extdata",package="rtop")
observations = readOGR(rpath,"observations")
predictionLocations = readOGR(rpath,"predictionLocations")

# Setting some parameters; nclus > 1 will start a cluster with nclus 
# workers for parallel processing
params = list(gDist = TRUE, cloud = FALSE, nclus = 1, rresol = 25)

# Create a column with the specific runoff:
observations$obs = observations$QSUMMER_OB/observations$AREASQKM

# Build an object
rtopObj = createRtopObject(observations, predictionLocations, 
           params = params)

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
