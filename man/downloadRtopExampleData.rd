\name{downloadRtopExampleData}
\alias{downloadRtopExampleData}
\title{
Download additional example data
}
\description{
Download additional example data from Vienna University of Technology
}
\usage{
downloadRtopExampleData(folder = system.file("extdata",package="rtop"))
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


\author{ Jon Olav Skoien }
\examples{
\dontrun{
  downloadRtopExampleData()
  rpath = system.file("extdata",package="rtop")
  setwd(rpath)
  observations = readOGR(".","observations") 
}
}
\keyword{plot}
