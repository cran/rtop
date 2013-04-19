\name{rtopDisc}
\alias{rtopDisc}
\alias{rtopDisc.rtop}
\alias{rtopDisc.SpatialPolygonsDataFrame}
\alias{rtopDisc.SpatialPolygons}
\alias{rtopDisc.rtopVariogram}
\title{Discretize areas}
\description{
\code{rtopDisc} will discretize an area for 
regularization or calculation of Ghosh-distance
}
\usage{
\method{rtopDisc}{rtop}(object, params = list(),...) 
\method{rtopDisc}{SpatialPolygonsDataFrame}(object, params = list(), bb = bbox(object), ...) 
\method{rtopDisc}{SpatialPolygons}(object, params = list(), bb = bbox(object),  ...) 
\method{rtopDisc}{rtopVariogram}(object, params = list(), ...) 
}
\arguments{
\item{object}{object of class \code{\link[sp]{SpatialPolygons}} or \code{\link[sp:SpatialPolygons]{SpatialPolygonsDataFrame}}
or \code{rtopVariogram},
or an object with class \code{rtop} that includes one of the above
}
\item{bb}{boundary box, usually modified to be the common boundary box for two 
          spatial object}
\item{params}{possibility to pass parameters to modify the standard parameters for 
      the \code{rtop} package, set in \code{\link{getRtopParams}}. Typical parameters
      to modify for this function are:
      \itemize{
        \item{rresol = 100; minimum number of discretization points in areas}  
        \item{hresol = 5; number of discretization points in one direction for areas in binned variograms}
        \item{hstype = "regular"; sampling type for binned variograms}
        \item{rstype = "rtop"; sampling type for real areas}
      }    
  }
\item{...}{Possibility to pass individual parameters}
} 

\value{ 
The function returns a list of discretized areas, or if called with an 
rtop-object as argument, the object with lists of discretizations of the
observations and prediction locations (if part of the object). If the function
is called with an rtopVariogram (usually this is an internal call), the 
list contains discretized pairs of hypothetical objects from each bin of 
the semivariogram with a centre-to-centre
distance equal to the average distance between the objects in a certain bin.
}

\details{
There are different options for discretizing the objects. When the areas
from the bins are discretized, the options are \code{random} or \code{regular} sampling, 
\code{regular} sampling is the default.

For the real areas, regular sampling appears to have computational advantages compared
with random sampling. In addition to the traditional regular sampling, \code{rtop}
also offers a third type of sampling which assures that the same discretization
points are used for overlapping areas. 

Starting with a coarse grid covering the 
region of interest, this will for a certain support be refined till a 
requested minimum number of points from the grid is within the support. 
In this way, for areal supports, the number of points in the area with the 
largest number of points will be approximately four times the requested minimum 
number of points. This methods also assure that points used to discretize a 
large support will be reused when discretizing smaller supports within the 
large one, e.g. subcatchments within larger catchments. 
}



\author{ Jon Olav Skoien }
\seealso{
\code{\link{rtop-package}} and \code{\link{rtopVariogram}}
}
\keyword{spatial}
