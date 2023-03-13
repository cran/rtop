\name{plot.rtopVariogramCloud}
\alias{plot.rtopVariogramCloud}
\title{
Plot and Identify Data Pairs on Sample Variogram Cloud
}
\description{
Plot a sample variogram cloud, possibly with identification of
individual point pairs
}
\usage{
\method{plot}{rtopVariogramCloud}(x, ...)
}
\arguments{
\item{x}{ object of class \code{variogramCloud}}
\item{...}{ parameters that are passed through to \code{\link[gstat]{plot.variogramCloud}}
The most important are:
\itemize{
\item{identify}{ logical; if TRUE, the plot allows identification of
a series of individual point pairs that correspond to individual variogram
cloud points (use left mouse button to select; right mouse button ends) }
\item{digitize}{ logical; if TRUE, select point pairs by digitizing a
region with the mouse (left mouse button adds a point, right mouse button
ends) }
\item{xlim}{ limits of x-axis }
\item{ylim}{ limits of y-axis }
\item{xlab}{ x axis label }
\item{ylab}{ y axis label }
\item{keep}{ logical; if TRUE and \code{identify} is TRUE, the labels
identified and their position are kept and glued to object x, which is
returned. Subsequent calls to plot this object will now have the labels
shown, e.g. to plot to hardcopy }
}

}
}
\note{
This function is mainly a wrapper around \code{\link[gstat]{plot.variogramCloud}},
necessary because of different column names and different class names. The 
description of arguments and value can therefore be found in
the help page of \code{\link[gstat]{plot.variogramCloud}}.}



\references{ 
\url{http://www.gstat.org/}
}
\author{ Jon Olav Skoien }
\seealso{
\link[gstat]{plot.gstatVariogram}
}
\examples{
\donttest{
rpath = system.file("extdata",package="rtop")
if (require("rgdal")) {
  observations = readOGR(rpath,"observations")
} else {
  library(sf)
  observations = st_read(rpath, "observations")
}
observations$obs = observations$QSUMMER_OB/observations$AREASQKM

# Create the sample variogram
rtopVario = rtopVariogram(observations, params = list(cloud = TRUE))
plot(rtopVario)

}
}

\keyword{dplot}
