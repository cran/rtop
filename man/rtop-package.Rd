\name{rtop-package}
\alias{rtop-package}

\title{A package providing methods for analysis and spatial interpolation
of data with an irregular support }

\description{
This package provides geostatistical methods for analysis and
interpolation of data that has an irregular support, such as as runoff
characteristics or population health data.
The methods in this package are based on the top-kriging approach
suggested in Skoien et al (2006), with some extensions from Gottschalk (1993).
This package can be used as an add-on package for the automatic interpolation
package developed within the intamap project (www.intamap.org).
}

\section{Workflow}{

The work flow within the package suggests that the user is interested in
a prediction of a process at a series of locations where observations have not been made.
The example below shows a regionalization of mean annual runoff in Austria.

Although it is possible to perform each step with all necessary arguments, 
the easiest interface to
the method is to store all variables (such as observations, prediction locations
and parameters) in an rtop-object, which is created by a call to
\code{\link{createRtopObject}}. The element \code{params} below consists of
changes to the default parameters. A further description can be found
in \code{\link{getRtopParams}}. The changes below means that
the functions will use geostatistical distance instead of full regularization,
and that the variogram model will be fitted to the variogram cloud.
Most other functions in the \code{rtop}-package can take this object as an argument,
and will add the results as one or more new element(s) to this object.

The data in the example below are stored as shape-files
in the extdata-directory of the rtop-pacakge, use the directory of your own data instead. 
The observations consist of mean summer runoff
from 138 catchments in Upper Austria. The predictionLocations are 863 catchments
in the same region. observations and predictionLocations are stored as
\code{\link[sp]{SpatialPolygonsDataFrame}}-objects.


\preformatted{
library(rgdal)
rpath = system.file("extdata",package="rtop")
observations = readOGR(rpath,"observations")

# Create a column with the specific runoff:
observations$obs = observations$QSUMMER_OB/observations$AREASQKM
predictionLocations = readOGR(rpath,"predictionLocations")
params = list(gDist = TRUE, cloud = TRUE)
rtopObj = createRtopObject(observations,predictionLocations, params = params)
}

There are help-methods available in cases when data are not available as
shape-files, or when the observations are not part of the shape-files.
See \code{\link{readAreaInfo}} and \code{\link{readAreas}}.

A call to \code{\link{rtopVariogram}} adds the sample variogram to the object,
whereas \code{\link{rtopFitVariogram}} fits a variogram model. The last function
will call \code{\link{rtopVariogram}} if \code{rtopObj} does not contain a sample variogram.

\preformatted{
rtopObj = rtopVariogram(rtopObj)
rtopObj = rtopFitVariogram(rtopObj)
}


The function \code{\link{checkVario}} is useful to produce some
diagnostic plots for the sample variogram and the fitted variogram model.

\preformatted{
checkVario(rtopObj)
}

The interpolation function (\code{\link{rtopKrige}}) solves the kriging system based on the
computed regularized semivariances. The covariance matrices are created in a
separate regularization function (\code{\link{varMat}}), and are stored in
the rtop-object for easier access if it is necessary to redo parts of the
analysis, as this is the computationally expensive part of the interpolation.
Cross-validation can be called with the argument \code{cv=TRUE}, either in
\code{params} or in the call to \code{\link{rtopKrige}}.


\preformatted{
rtopObj = rtopKrige(rtopObj)
spplot(rtopObj$predictions,col.regions = bpy.colors(), c("var1.pred","var1.var"))
rtopObj = rtopKrige(rtopObj, cv = TRUE)
spplot(rtopObj$predictions,col.regions = bpy.colors(), c("var1.pred","var1.var"))
}



}
\references{
L. Gottschalk. Interpolation of runoff applying objective methods.
Stochastic Hydrology and Hydraulics, 7:269-281, 1993.

Skoien J. O., R. Merz, and G. Bloschl. Top-kriging - geostatistics on stream networks.
Hydrology and Earth System Sciences, 10:277-287, 2006.
}
\keyword{spatial}

