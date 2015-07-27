\name{readAreas}
\alias{readAreas}
\title{ help file for creating SpatialPolygonsDataFrame with observations and/or predictionLocations
of data with a spatial support }
\description{
readAreas will read area-files, add observations and convert the result to \cr
\code{\link[sp]{SpatialPolygonsDataFrame}}
}
\usage{
readAreas(object, adir=".",ftype = "xy",projection = NA, ...) 
}
\arguments{
\item{object}{either name of file with areal information or \code{\link[sp]{SpatialPointsDataFrame}}
   with observations}
\item{adir}{directory where the files with areal information are to be found}
\item{ftype}{type of file, the only type supported currently is "xy", 
   referring to x- and y-coordinates of boundaries}
\item{projection}{add projection to the object if input is boundary-files }
\item{...}{further parameters to be passed to \code{\link{readAreaInfo}}}
} 

\value{ 
The function creates a \code{\link[sp:SpatialPolygons]{SpatialPolygonsDataFrame}} of observations
and/or predictionLocations, depending on the information given in \code{object}.
}


\details{
If \code{object} is a file name, \code{\link{readAreaInfo}} will be called.
If it is a \cr
\code{\link[sp:SpatialPoints]{SpatialPointsDataFrame}} with observations and/or 
predictionLocations, the function will read areal data from files according 
to the ID associated with each observation/predictionLocation. 

The function is of particular use when data are not available as
shape-files, or when the observations are not part of the shape-files.
This function is mainly for compatibility with the former FORTRAN-version.
The simplest way to read the data in that case is through \code{readShapePoly} 
in the \code{maptools}-package or \code{readOGR} in 
the \code{rgdal}-package. See
also \code{\link{rtop-package}}.
}

\author{ Jon Olav Skoien }
\keyword{spatial}
