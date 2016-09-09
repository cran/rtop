\name{getRtopParams}
\alias{getRtopParams}
\title{ Setting parameters for the intamap package}
\description{ This function sets a range of the parameters for the intamap package, 
to be included in the object described in \code{\link{rtop-package}}}

\usage{getRtopParams(params,newPar, observations, formulaString, ...)
}
\arguments{
\item{params}{An existing set of parameters for the interpolation process, 
              of class \cr
              \code{intamapParams} or a list of parameters for modification
              of the default parameters}
\item{newPar}{A \code{list} of parameters for updating \code{params} or for
              modification of the default parameters. 
              Possible parameters with their defaults are given below}
\item{observations}{\code{\link[sp]{SpatialPolygonsDataFrame}} with observations,
                 used for setting some of the default parameters}
\item{formulaString}{formula that defines the dependent variable as a linear model 
      of independent variables, see e.g. \code{\link{createRtopObject}} for more details.}
\item{...}{Individual parameters for updating \code{params} or for
              modification of the default parameters.
              Possible parameters with their defaults are given below 



  \itemize{
     \item{model = "Ex1"}{ - variogram model type. Currently the following models are implemented:
          \itemize{
              \item{Exp}{ - Exponential model}
              \item{Ex1}{ - Multiplication of a modified exponential and fractal model, 
                          the same model as used in Skoien et al(2006).}
              \item{Gau}{ - Gaussian model}
              \item{Ga1}{ - Multiplication of gaussian and fractal model}
              \item{Sph}{ - Spherical model}
              \item{Sp1}{ - Multiplication of spherical and fractal model}
              \item{Fra}{ - Fractal model}
          }}
     \item{parInit}{ - the initial parameters and the limits of the variogram model to be fitted,
           given as a matrix with three columns, where the first column is the 
           lower limit, the second column is the upper limit and the third column
           are starting values.}
     \item{nugget = FALSE}{ - logical; if TRUE, nugget effect should be estimated}
     \item{unc = TRUE}{ - logical; if TRUE the variance of observations are in column \code{unc}}
     \item{rresol = 100}{ - minimum number of discretization points in each area}  
     \item{hresol = 5}{ - number of discretization points in one direction for elements in binned variograms}
     \item{cloud = FALSE}{ - logical; if TRUE use the cloud variogram for variogram fitting}
     \item{amul = 1}{ - defines the number of areal bins within one order of magnitude. Numbers between 1 and 3 
                      are possible, as this parameter refers to the \code{axp} parameter of
                      \code{\link[graphics]{axTicks}}.}
     \item{dmul = 3}{ - defines the number of distance bins within one order of magnitude. Numbers between 1 and 3 
                      are possible, as this parameter refers to the \code{axp} parameter of
                      \code{\link[graphics]{axTicks}}.}
     \item{fit.method = 8}{ - defines the type of Least Square method for fitting of variogram. 
                  The methods 1-7 correspond to the similar methods in \code{\link{gstat}}.
       \itemize{
        \item{1}{ - weighted least squares with number of pairs per bin: \cr 
                    err = n * (yobs-ymod)^2}
        \item{2}{ - weighted least squares difference according to Cressie (1985): \cr
                    err2 = abs(yobs/ymod-1)}
        \item{6}{ - ordinary least squares difference: err = (yobs-ymod)^2}
        \item{7}{ - similar to default of gstat, where higher weights are given to shorter distances err = n/h^2 * (yobs-mod)^2}
        \item{8}{ - Opposite of weighted least squares difference according to Cressie (1985): err3=abs(ymod/yobs-1)}
        \item{9}{ - neutral WLS-method - err = min(err2,err3)}
       }
     }
     \item{gDistEst = FALSE}{ - use geostatistical distance when fitting variograms}
     \item{gDistPred = FALSE}{ - use geostatistical distance for semivariogram matrices}
     \item{gDist}{ - parameter to set jointly \code{gDistEst = gDistPred = gDist}}
     \item{nmax = 10}{for local kriging: the number of nearest observations that
            should be used for a kriging prediction or simulation, where
            nearest is defined in terms of the space of the spatial locations.
            By default, 10 observations are used.}
     \item{maxdist = Inf}{ - for local kriging: only observations within a distance 
                      of \code{maxdist} from the prediction location are used for prediction 
                      or simulation; if combined with nmax, both criteria apply }
     \item{hstype = "regular"}{ - sampling type for binned variograms}
     \item{rstype = "rtop"}{ - sampling type for the elements, see also \code{\link{rtopDisc}}}
     \item{nclus = 1}{- number of CPUs to use if parallel processing is wanted; nclus = 1 means
                       no parallelization}
     \item{cnAreas = 100}{- limit whether parallel processing should be applied; the minimum number
                       of areas in \code{\link{varMat}}, and also controlling when to use parallel 
                       processing in
                       \code{\link{rtopDisc}}, when \cr
                       \code{nAreas*params$rresol/100 > cnAreas} }
     \item{clusType = NULL}{- the cluster type to be started for parallel processing; uses the default 
                       type of the system when clusType = NULL}
     \item{outfile = NULL}{file where output can be printed during parallel execution}
     \item{varClean = FALSE}{logical; if TRUE it will remove highly correlated areas from the covariance matrix during simulation }
     \item{wlim = 1.5}{ - an upper limit for the norm of the weights in kriging, see \code{\link{rtopKrige}}}
     \item{wlimMethod = "all"}{which method to use for reducing the norm of the weights if necessary. Either "all", which modifies all weights equally or "neg" which reduces negative weights and large weights more than the smallest weights }
     \item{cv = FALSE}{ - logical; for cross-validation of observations}
     \item{debug.level = 1}{ - used in some functions for giving additional output. See 
                individual functions for more information.}
     \item{partialOverlap = FALSE}{whether to work with partially overlapping areas}
     \item{olim = 1e-4}{smallest overlapping area to be used for partial overlap, relative to the smallest of the areas}
     \item{nclus = 1}{option to use parallel processing, nclus > 1 defines the number of workers to be started}
     \item{clusType = NA}{which cluster type to start if nclus > 1; the default is used if nclusType = NA }
     \item{cnAreas = 200}{The minimum number of observations or observations plus predictions allowing parallelization in
                    the creation of the covariance matrix}
     \item{cDlim = 1e6}{The minimum number of discretization points for allowing parallelization in the discretization process}
     \item{observations}{ - used for initial values of parameters if supplied}
     \item{formulaString}{ - used for initial values of parameters if supplied}
}
}
}
\value{ 
A list of the parameters with class \code{rtopParams} to be included in the 
\code{object} described in \link{rtop-package}}




\note{
This function will mainly be called by \code{\link{createRtopObject}}, but 
can also be called by the user to create a parameter set or update an 
existing parameter set. If none of the arguments is a list of class
\code{rtopParams}, the function will assume that the argument(s) are
modifications to the default set of parameters. The function can also be called
by other functions in the rtop-package if the users chooses not to work with
an object of class \code{rtop}.

If the function is called with two lists of parameters (but the first one is
not of class \code{rtopParams}) they are both seen as modifications to the 
default parameter set. If they share some parameters, the parameter values from
the second list will be applied.

Parallel processing has been included for some of the functions. The default is no 
parallel procesing, and the package also attempts to decide whether it is sensible to 
start a set of clusters and distribute jobs to them based on the size of the job. 
The default limit might not be the best for every system.


}
\references{ 
Cressie, N. 1985. Fitting variogram models by weighted least squares. Mathematical Geology, 17 (5), 563-586

Skoien J. O., R. Merz, and G. Bloschl. Top-kriging - geostatistics on stream networks.
Hydrology and Earth System Sciences, 10:277-287, 2006

Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione, A., 2014. Rtop: An R package for interpolation of data with a variable spatial support, with an example from river networks. Computers & Geosciences, 67.
}
\author{ Jon Olav Skoien }

\seealso{
\code{\link{createRtopObject}} and \code{\link{rtop-package}}
}
\examples{
# Create a new set of intamapParameters, with default parameters:
params = getRtopParams()
# Make modifications to the default list of parameters
params = getRtopParams(newPar = list(gDist = TRUE, nugget = FALSE))
# Make modifications to an existing list of parameters
params = getRtopParams(params = params, newPar = list(gDist = TRUE,
         nugget = FALSE))
}
\keyword{spatial}
