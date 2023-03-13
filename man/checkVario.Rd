\name{checkVario}
\alias{checkVario}
\alias{checkVario.rtop}
\alias{checkVario.rtopVariogramModel}
\title{ Plot variogram fitted to data with support }
\description{
The function will create diagnostic plots for analysis of the variograms fitted
to sample variograms of data with support
}
\usage{
\method{checkVario}{rtop}(object, acor = 1, log = "xy", cloud = FALSE, 
           gDist = TRUE, acomp = NULL, curveSmooth = FALSE, params = list(), ...) 

\method{checkVario}{rtopVariogramModel}(object, 
           sampleVariogram = NULL, observations = NULL, 
           areas = NULL, dists = NULL, acomp = NULL,  
           params = list(), compVars = list(), acor = 1, 
           log = "xy", legx = NULL, legy = NULL, 
           plotNugg = TRUE, curveSmooth = FALSE, ...) 
}



\arguments{
\item{object}{either: object of class \code{rtop} (see \code{\link{rtop-package}}), or an object of type \cr
\code{rtopVariogram} }
\item{acor}{unit correction factor in the key, e.g. to see numbers more easily 
        interpretable for large areas. As an example, ucor = 0.000001 when area
        is given in square meters and should rather be shown as square kilometers.
        Note that this parameter also changes the value of the nugget to the new unit.}
\item{log}{text variable for log-plots, default to log-log \code{"xy"}, can otherwise be set
        to \code{"x"}, \code{"y"} or \code{""}}
\item{cloud}{logical; whether to look at the cloud variogram instead of the binned variogram}
\item{gDist}{logical; whether to use ghosh-distance for semivariogram regularization instead of full integration of the semivariogram}
\item{sampleVariogram}{a sample variogram of the data}
\item{observations}{a set of observations}
\item{areas}{either an array of areas that should be used as examples, or 
        the number of areas per order of magnitude (similar to the parameter \code{amul};
         see \code{\link{getRtopParams}}. amul from \code{rtopObj}
        or from the default parameter set will be used if not defined here.}
\item{dists}{either an array of distances that should be used as examples, or
        the number of distances per order of magnitude(similar to the parameter \code{amul};
        see \code{\link{getRtopParams}}. amul from \code{rtopObj}
        or from the default parameter set will be used if not defined here.}
\item{acomp}{either a matrix with the area bins that should be visualized, or a number 
        giving the number of pairs to show. If a sample variogram is given, the \code{acomp}
        pairs with highest number of pairs will be used}
\item{curveSmooth}{logical or numerical; describing whether the curves in the last plot should be smoothed or not. If numeric,
        it gives the degrees of freedom (df) for the splines used for smoothing. See also \code{\link[stats]{smooth.spline}}}
\item{params}{list of parameters to modify the default parameters of rtopObj or 
        the default parameters found from \code{\link{getRtopParams}}}
\item{compVars}{a list of variograms of \code{gstat}-type for comparison, see 
        \code{\link[gstat:vgm]{vgm}}. The names of the variograms in the list will
        be used in the key.}

\item{legx}{x-coordinate of the legend for fine-tuning of position, see x-argument of \cr
\code{\link[graphics]{legend}}}
\item{legy}{y-coordinate of the legend for fine-tuning of position, see y-argument of \cr
\code{\link[graphics]{legend}}}
\item{plotNugg}{logical; whether the nugget effect should be added to the plot or not}
\item{...}{arguments to lower level functions}

} 

\value{ 
The function gives diagnostic plots for the fitted variograms, where the
regularized variograms are shown together with the sample variograms and 
possibly also user defined variograms. In addition, if an rtopObject
is submitted, the function will also give plots of the relationship between
variance and area size and a scatter plot of the fit of the observed and
regularized variogram values. The sizes of the dots are relative to the number
of pairs in each group.
}

\references{ 
Skoien J. O., R. Merz, and G. Bloschl. Top-kriging - geostatistics on stream networks.
Hydrology and Earth System Sciences, 10:277-287, 2006.

Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione, A., 2014. Rtop: An R package for interpolation of data with a variable spatial support, with an example from river networks. Computers & Geosciences, 67.
}


\author{ Jon Olav Skoien }
\examples{
\donttest{
library(gstat)
rpath = system.file("extdata",package="rtop")
if (require("rgdal")) {
  observations = readOGR(rpath,"observations")
  predictionLocations = readOGR(rpath,"predictionLocations")
} else {
  library(sf)
  observations = st_read(rpath, "observations")
  predictionLocations = st_read(rpath,"predictionLocations")
}

# Create a column with the specific runoff:
observations$obs = observations$QSUMMER_OB/observations$AREASQKM
params = list(cloud = TRUE, gDist = TRUE)
rtopObj = createRtopObject(observations, predictionLocations, 
                           params = params)

# Fit a variogram (function also creates it)
rtopObj = rtopFitVariogram(rtopObj)
checkVario(rtopObj, 
    compVar = list(first = vgm(5e-6, "Sph", 30000,5e-8), 
                   second = vgm(2e-6, "Sph", 30000,5e-8)))

rtopObj = checkVario(rtopObj, acor = 0.000001, 
          acomp = data.frame(acl1 = c(2,2,2,2,3,3,3,4,4), 
          acl2 = c(2,3,4,5,3,4,5,4,5)))
rtopObj = checkVario(rtopObj, cloud = TRUE, identify = TRUE, 
          acor = 0.000001)
}
}
\keyword{spatial}
