\name{sceua}
\alias{sceua}
\title{
Optimisation with the Shuffle Complex Evolution method}
\description{
Function to minimize the value of an objective function for calibration
}
\usage{
sceua(OFUN, pars, lower, upper, maxn = 10000, kstop = 5, pcento = 0.01,
    ngs = 5, npg = 5, nps = 5, nspl = 5, mings = 5, iniflg = 1, iprint = 0, iround = 3, 
    peps = 0.0001, plog = rep(FALSE,length(pars)), implicit = NULL, ...)
}
\arguments{
\item{OFUN}{ A function to be minimized, with first argument the vector of 
          parameters over which minimization is to take place. It should return a scalar result 
          as an indicator of the error for a certain parameter set}
\item{pars}{ a vector with the initial guess the parameters }
\item{lower}{ the lower boundary for the parameters }
\item{upper}{ the upper boundary for the parameters }
\item{maxn}{ the maximum number of function evaluations}
\item{kstop}{number of shuffling loops in which the criterion value must change by the given percentage before optimization is terminated}
\item{pcento}{percentage by which the criterion value must change in given number (kstop) of shuffling loops to continue optimization}
\item{ngs}{number of complexes in the initial population}
\item{npg}{number of points in each complex}
\item{nps}{number of points in a sub-complex}
\item{nspl}{number of evolution steps allowed for each complex before complex shuffling}
\item{mings}{minimum number of complexes required, if the number of complexes is allowed to reduce as the optimization proceeds}
\item{iniflg}{flag on whether to include the initial point in population
 = 0, not included
 = 1, included
}
\item{iprint}{flag for controlling print-out after each shuffling loop
     = 0, print information on the best point of the population
      = 1, print information on every point of the population
}
\item{iround}{number of significant digits in print-out}
\item{peps}{convergence level for parameter set (lower number means smaller difference between parameters of the population required for stop)}
\item{plog}{whether optimization should be done in log10-domain. Either a single TRUE value for all parameters, or a vector with TRUE/FALSE for the different parameters}
\item{implicit}{function for implicit boundaries for the parameters (e.g. sum(pars[4]+pars[5]) < 1). See below for details}
\item{...}{arguments for the objective function, must be named}
}







\value{
The function returns a list with the following elements
\itemize{
  \item{par}{a vector of the best parameters combination }
  \item{value}{value of the objective function for this parameter set}
  \item{convergence}{list of two values \itemize{
       \item{funConvergence}{function convergence relative to pcento}
       \item{parConvergence}{parameter convergence relative to peps}
       }}
  \item{counts}{number of function evaluations}
  \item{iterations}{number of shuffling loops}
}
}

\details{
sceua is an R-implementation of the Shuffle Complex Evolution - University of Arizona (Duan et al., 1992), 
a global optimization method which "combines the strengths of the simplex procedure of Nelder and Mead (1965)
with the concepts of controlled random search (Price, 1987), competetive evolusion (Holland, 1975)" with the 
concept of complex shuffling, developed by Duan et al. (1992). 


This implementation follows the Fortran implementation relatively close, but adds the possibility of searching 
in log-space for one or more of the parameters, and it uses the capability of R to pass functions
as arguments, making it possible to pass implicit conditions to the parameter selection. 

The objective function \code{OFUN} is a function which should give an error value for each parameter set.
It should never return non-numeric values such as NA, NULL, or Inf. If some parameter combinations can 
give such values, the return value should rather be a large number.

The function works with fixed upper and lower boundaries for the parameters. If the possible range of 
a parameter might span several orders of magnitude, it might be better to search in log-space for the optimal parameter, 
to reduce the risk of being trapped in local optima. This can be set with the argument \code{plog}, which is either
a single value (FALSE/TRUE) or a vector for all parameters. 
\code{plog = c(TRUE, FALSE, FALSE, TRUE, TRUE)} means that the search for parameters 1,4 and 5 should be in log10-space, 
whereas the search for parameters 2 and 3 are in normal space.

Implicit boundaries can be evoked by passing a function \code{implicit} to \code{sceua}.
This function should give 0 when parameters are acceptable 
and 1 if not. If, for example, the condition is that the following sum of parameters four and five should be limited:

sum(pars[4]+pars[5]) <= 1

then the function will be implicit = function(pars) (2*pars[4] + pars[5]) > 1
}

\references{Duan, Q., Sorooshian, S., and Gupta, V.K., 1992. Effective and efficient 
global optimization for conceptual rainfall-runoff models.
Water Resour. Res. 28 (4), 1015?1031.

Holland, H.H., 1975. Adaptation in natural and artificial systems, 
University of Michigan Press, Ann Arbor.

Nelder, J.A. and Mead, R., 1965. A simplex method for function minimization, 
Comput. J., 7(4), 308-313.

Price, W.L., 1987. Global optimization algorithms for a CAD workstation, 
J. Optim. Theory Appl., 55(1), 133-146.



}

\author{ Jon Olav Skoien }
\examples{
set.seed(1)
# generate example data from a function with three parameters
# with some random noise
fun = function(x, pars) pars[2]*sin(x*pars[1])+pars[3]
x = rnorm(50, sd = 3)
y = fun(x, pars = c(5, 2, 3)) +  rnorm(length(x), sd = 0.3)
plot(x,y)
   
# Objective function, summing up squared differences
OFUN = function(pars, x, yobs) {
  yvals = fun(x, pars)  
  sum((yvals-yobs)^2)
}

sceuares = sceua(OFUN, pars = c(0.1,0.1,0.1), lower = c(-10,0,-10), 
                 upper = c(10,10,10), x = x, yobs = y)
sceuares
xx = seq(min(x), max(x), 0.1)
lines(xx, fun(xx, pars = sceuares$par))

}