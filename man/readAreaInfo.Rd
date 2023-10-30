\name{readAreaInfo}
\alias{readAreaInfo}
\title{ create SpatialPointsDataFrame with observations
of data with a spatial support }
\description{
readAreaInfo will read a text file with observations and descriptions of 
data with a spatial support.
}
\usage{
readAreaInfo(fname = "ainfo.txt", id = "id",
                iobs = "iobs", obs = "obs", unc = "unc",
                filenames= "filenames", sep = "\t", 
                debug.level = 1, moreCols = list(NULL)) 
}
\arguments{
\item{fname}{ name of file with areal information}
\item{id}{ name of column with observation id }
\item{iobs}{ name of column with number of observations}
\item{obs}{ name of column with observations }
\item{unc}{ name of column with possible uncertainty of observation}
\item{filenames}{ name of column with filenames of areas if different names than id should
be used.}
\item{sep}{separator in csv-file}
\item{debug.level}{used for giving additional output}
\item{moreCols}{name of other column names the user wants included in ainfo}
} 

\value{ 
\code{\link[sp:SpatialPoints]{SpatialPointDataFrame}} with information about observations and/or
predictionLocations.
}

\details{
The function is of particular use when data are not available as
shape-files, or when the observations are not part of the shape-files.
This function is mainly for compatibility with the former FORTRAN-version.
The simplest way to read the data in that case is through \code{\link[sf]{st_read}}. See
also \code{\link{rtop-package}}.
}

\author{ Jon Olav Skoien }
\keyword{spatial}
