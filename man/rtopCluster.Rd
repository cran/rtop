\name{rtopCluster}
\alias{rtopCluster}
\title{start, access, stop or restart a cluster for parallel computation with rtop}
\description{Convenience function for using parallel computation with rtop. The function is usually not 
called by the user.
}
\usage{
rtopCluster(nclus, ..., action = "start", type, outfile = NULL)
}

\arguments{
\item{nclus}{The number of workers in the cluster}
\item{...}{Arguments for \code{\link[parallel]{clusterEvalQ}}; commands to be evaluated for each worker, such as library-calls}
\item{action}{Defines the action of the function. There are three options:
\describe{
 \item{"start"}{Starts a new cluster if necessary, reuses an existing if it has already been started}
 \item{"restart"}{Stops the cluster and starts it again. To be used in case there are difficulties with the cluster, 
     or if the user wants to change the type of the cluster}
}
}
\item{type}{The type of cluster; see \code{\link[parallel]{makeCluster}} for more details. 
            The default of makeCluster is used if type is missing or NA.}
\item{outfile}{File to direct the output, \code{\link[parallel]{makeCluster}} for more details. }
}
 

\value{ 
If the function is called with action equal to "start" or "restart", the result is a cluster with nclus workers.
The cluster is also added to the global options with the name rtopCluster \cr
(\code{getOption("rtopCluster")}).

If the function is called with action equal to "stop", the function stops the cluster, sets the rtopCluster of options to NULL and returns NULL to the user.


}

\details{
It is usually not necessary for the user to call this function for starting or accessing a cluster. This is done automatically
by the different rtop-functions when needed if the parameter nclus is larger than one (see \code{\link{getRtopParams}}). If the user actually starts the cluster by a call to this function, it will 
also be necessary to set the nclus parameter to a value larger than one for the cluster to be used by different functions.

Restarting the cluster might be necessary if the cluster has a problem (e.g. does not return memory) or if the user wants to change to a different cluster type. 

Stopping the cluster is useful when the user does not want to continue with parallel computation and wants to close down the workers.

}



\author{ Jon Olav Skoien }

\keyword{spatial}
