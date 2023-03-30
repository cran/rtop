\name{rtopSim}
\alias{rtopSim}
\alias{rtopSim.rtop}
\alias{rtopSim.default}
\title{Spatial simulation of data with spatial support}
\description{
rtopSim will conditionally or unconditionally simulate data with areal support. This function should be seen as experimental, some issues are described below.
}
\usage{


  
\method{rtopSim}{rtop}(object, varMatUpdate = FALSE, beta = NA, 
      largeFirst = TRUE, replace = FALSE, params = list(), 
      dump = NULL, debug.level, ...) 
\method{rtopSim}{default}(object = NULL, predictionLocations,
    varMatObs,  varMatPredObs, varMatPred, variogramModel, ...) 
}

\arguments{
\item{object}{object of class \code{rtop} or \code{\link[sp]{SpatialPolygonsDataFrame}} or
             \code{sf} (\code{\link[sf]{st_sf}}) or \code{NULL}}
\item{varMatUpdate}{logical; if TRUE, also existing variance matrices will 
            be recomputed, if FALSE, only missing variance matrices will be computed,
            see also \code{\link{varMat}}}
\item{beta}{The expected mean of the data, for unconditional simulations}
\item{largeFirst}{Although the simulation method follows a random path around the 
            predictionLocations, simulating the largest area first will assure
            that the true mean of the simulated values will be closer to beta}
\item{replace}{logical; if observation locations are also present as predictions, should they be replaced?
               This is particularly when doing conditional simulations for a set of observations 
               with uncertainty.}
\item{params}{a set of parameters, used to modify the standard parameters for 
      the \code{rtop} package, set in \code{\link{getRtopParams}}. The argument params can also be used
      for the other methods, through the ...-argument.}
\item{dump}{file name for saving the updated object, after adding variance matrices. Useful if there
      are problems with the simulation, particularly if it for some reason crashes.}
\item{debug.level}{debug.level} {logical; controls some output, will override the object parameters}
\item{predictionLocations}{\code{\link[sp]{SpatialPolygons}} or \code{\link[sp]{SpatialPolygonsDataFrame}}
     or \code{sf}-object with locations for simulations. }
\item{varMatObs}{covariance matrix of possible observations, where diagonal must consist 
      of internal variance, typically generated from call
      to \code{\link{varMat}} }
\item{varMatPredObs}{covariance matrix between possible observation locations and 
      simulation locations, typically generated from call
      to \code{\link{varMat}}}
\item{varMatPred}{covariance matrix between simulation locations, typically generated
      from a call to \code{\link{varMat}}}
\item{variogramModel}{a variogram model of type \code{\link{rtopVariogramModel}} }
\item{...}{possible modification of the object parameters or default parameters.}
}

\value{ 
If called with \code{\link[sp]{SpatialPolygons}} or \code{sf} as predictionLocations
and either \cr
\code{\link[sp]{SpatialPolygonsDataFrame}}, \code{sf} or \code{NULL} for
observations, the function returns a\cr
\code{\link[sp]{SpatialPolygonsDataFrame}} or \code{sf} with simulations at the 
locations defined in \cr
\code{predictionLocations}

If called with an rtop-object, the function returns the same object with the
simulations added to the object.
}

\details{
This function can do constrained or unconstrained simulation for areas.
The simplest way of calling the function is with an rtop-object that 
contains the fitted variogram model and all the other necessary data (see
\code{\link{createRtopObject}} or \code{\link{rtop-package}}). \code{rtopSim}
is the only function in \code{rtop} which does not need  observations. 
However, a variogram model is still necessary to perform simulations.

The arguments \code{beta} and \code{largeFirst} are only used for unconditional simulations.


The function is still in an experimental stage, and might change in the future. There are some issues with the current implementation:
\itemize{
\item{Numerical issues can in some cases give negative estimation variances, which will result in an invalid distribution for the simulation. This will result in simulated NA values for these locations.}
\item{The variability of simulated values for small areas (such as small headwater catchments) will be relatively high based on the statistical uncertainty. This could be overestimated compared to the uncertainty which is possible based on rainfall.}
}
}




\references{ 
Skoien J. O., R. Merz, and G. Bloschl. Top-kriging - geostatistics on stream networks.
Hydrology and Earth System Sciences, 10:277-287, 2006.

Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione, A., 2014. Rtop: An R package for interpolation of data with a variable spatial support, with an example from river networks. Computers & Geosciences, 67.

}
\author{ Jon Olav Skoien }

\examples{
\donttest{
# The following command will download  the complete example data set
# downloadRtopExampleData() 

rpath = system.file("extdata",package="rtop")
if (require("rgdal") & FALSE) { # Preparing for retirement of rgdal
  observations = readOGR(rpath,"observations")
  predictionLocations = readOGR(rpath,"predictionLocations")
} else {
  library(sf)
  observations = st_read(rpath, "observations")
  predictionLocations = st_read(rpath,"predictionLocations")
}

# Setting some parameters; nclus > 1 will start a cluster with nclus 
# workers for parallel processing
params = list(gDist = TRUE, cloud = FALSE, nclus = 1, rresol = 25)

# Create a column with the specific runoff:
observations$obs = observations$QSUMMER_OB/observations$AREASQKM

# Build an object
rtopObj = createRtopObject(observations, predictionLocations, 
                           params = params, formulaString = "obs~1")

# Fit a variogram (function also creates it)
rtopObj = rtopFitVariogram(rtopObj)

# Conditional simulations for two new locations 
rtopObj10 = rtopSim(rtopObj, nsim = 5)
rtopObj11 = rtopObj

# Unconditional simulation at the observation locations
# (These are moved to the predictionLocations)
rtopObj11$predictionLocations = rtopObj11$observations
rtopObj11$observations = NULL
# Setting varMatUpdate to TRUE, to make sure that covariance
# matrices are recomputed
rtopObj12 = rtopSim(rtopObj11, nsim = 10, beta = 0.01, 
                    varMatUpdate = TRUE)

summary(data.frame(rtopObj10$simulations))
summary(data.frame(rtopObj12$simulations))

}
}
\keyword{spatial}
