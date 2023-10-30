\name{variogramModel}
\alias{rtopVariogramModel}
\alias{updateRtopVariogram}
\alias{updateRtopVariogram.rtop}
\alias{updateRtopVariogram.rtopVariogramModel}
\title{ create or update variogram model}
\description{
This gives an easier interface to the parameters of the variogram model
}
\usage{
rtopVariogramModel(model = "Ex1", sill = NULL, range = NULL, 
    exp = NULL, nugget = NULL, exp0 = NULL,
    observations = NULL, formulaString = obs~1) 
\method{updateRtopVariogram}{rtop}(object, ...)
\method{updateRtopVariogram}{rtopVariogramModel}(object, action = "mult", ..., 
              checkVario = FALSE, 
sampleVariogram = NULL, observations = NULL) 

}


\arguments{
\item{model}{variogram model, currently "Ex1" is the only implemented, see Skoien et al (2006)}
\item{sill}{sill of variogram}
\item{range}{range of variogram}
\item{exp}{the exponent of the fractal part of the variogram, see Skoien et al (2006)}
\item{exp0}{gives the angle of the first part of the variogram in a log-log plot (weibull type), 
      should be between 0 and 2. See Skoien et al (2006)}
\item{nugget}{nugget of point variogram}
\item{formulaString}{formula that defines the dependent variable as a linear model 
      of independent variables, see e.g. \code{\link{createRtopObject}} for more details.}
\item{object}{either: object of class \code{rtop} (see \link{rtop-package}), or an rtopVariogramModel. }
\item{action}{character variable defining whether the new parameters should 
      be \code{add}(-ed), \code{mult}(-iplied) or \code{replace} the former parameters.
      Leaving the parameters equal to NULL will cause no change.}
\item{checkVario}{logical, will issue a call to\code{\link{checkVario}} if TRUE}
\item{sampleVariogram}{a sample variogram of the data}
\item{observations}{a set of observations}
\item{...}{parameters to lower level functions}

} 

\value{ 
The function helps creating and updating the parameters of the variogram, by 
using common names and simple update methods. This is mainly for manual 
fitting of the variogram. The automatic call to checkVario makes it easier to 
visualize the effect of the changes to the variogram}


\author{ Jon Olav Skoien }
\seealso{\code{\link{rtop-package}}}
\examples{
\dontrun{
library(sf)
rpath = system.file("extdata",package="rtop")
observations = st_read(rpath,"observations")
# Create a column with the specific runoff:
observations$obs = observations$QSUMMER_OB/observations$AREASQKM
predictionLocations = st_read(rpath,"predictionLocations")
rtopObj = createRtopObject(observations,predictionLocations)
 # Fit a variogram (function also creates it)
rtopObj = rtopFitVariogram(rtopObj)
rtopObj = updateRtopVariogram(rtopObj, exp = 1.5, action = "mult", 
              checkVario = TRUE)
}
}
\keyword{spatial}
