\name{downloadRtopExampleData}
\alias{downloadRtopExampleData}
\title{
Download additional example data
}
\description{
Download additional example data from Vienna University of Technology
}
\usage{
downloadRtopExampleData(folder = system.file("extdata",
                        package="rtop"))
}
\arguments{
\item{folder}{the folder to which the downloaded data set will be copied}
}
\value{
The function will have as a side effect that additional example data is
downloaded from Vienna University of Techology. This will for the default
case replace the existing example data-set in the \code{rtop} package. Alternatively
the user can specify a separate directory for the data set.
}

\references{ 
Skoien J. O., R. Merz, and G. Bloschl. Top-kriging - geostatistics on stream networks.
Hydrology and Earth System Sciences, 10:277-287, 2006.

Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione, A., 2014. Rtop: An R package for interpolation of data with a variable spatial support, with an example from river networks. Computers & Geosciences, 67.
}


\author{ Jon Olav Skoien }
\examples{
\dontrun{
  downloadRtopExampleData()
  rpath = system.file("extdata",package="rtop")
  observations = readOGR(rpath,"observations") 
}
}
\keyword{plot}
