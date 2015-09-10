\name{iHaz}
\alias{iHaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of montone hazard and survival functions
}
\description{
The main function for estimating the hazard function (and consequently the survival function) for left truncated interval censored data under the assumption of a monotone hazard function. The inputs for this function are baseline times (left truncation time) and the two time points for interval censoring. The censoring interval can also be \code{(0,b)} or \code{(a,Inf)} to denote left and right censoring, respectively. Similarly, setting the truncation time as zero allows us to specialize to the case of censored data without left truncation.

}
\usage{
iHaz(data, ini.index = 1:3, inter.maxiter = 1000, inter.tol = 1e-04,
    main.maxiter = 1000, main.tol = 0.001,
    check.condition = c("derv", "KKT"), verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame or list with the following components: \code{t} basleine/left truncation times, \code{a} and \code{b} the endpoints for interval censoring with \code{a < x < b} where \code{x} is the true (unobserved) event time.
}
  \item{ini.index}{
The initial index set for the support reudction algorithm. See details below.
}
  \item{inter.maxiter}{
The maximum number of iterations for secondary/intermediate algorithm for optimization over the reduced support set. See details below.
}
  \item{inter.tol}{
The tolerance to determine the stopping condition for the intermediate algorithm. See details below.
}
  \item{main.maxiter}{
The maximum number of iterations for the main loop of the algorithm. See details below.
}
  \item{main.tol}{
The tolerance to determine the stopping condition for the main loop. See details below.
}
  \item{check.condition}{
A string to specify the type of convergence criteria. The possible options are "KKT" which uses the Karush
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}
\details{
The estimated hazard function via the projection algorithm of Pan et al. (1998) is a step function with a potential change point at every time point in our dataset. In practice however, the estimated hazard has only a small number of change points.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, ini.index = 1:3, inter.maxiter = 1000, inter.tol = 1e-04,
    main.maxiter = 1000, main.tol = 0.001, check.condition = c("derv",
        "KKT"), verbose = FALSE)
{
    index <- ini.index
    if (verbose)
        cat("Building D matrix...")
    a <- data$a
    t <- data$t
    b <- data$b
    z <- sort(unique(c(a, b, t)))
    Dmatrix <- cpp_buildMatrixD(t, a, b, z)
    if (verbose)
        cat("done.\n\n")
    for (i in 1:main.maxiter) {
        if (verbose)
            cat(paste0("Iteration Number: ", i, "\n"))
        myans <- InterAlg(ini = 1, data, index = index, epsilon = inter.tol,
            maxiter = inter.maxiter, Dmatrix = Dmatrix)
        if (check.condition[1] == "derv") {
            temp <- check.derv(myans)
            newindex <- cpp_findLocalMax(temp, main.tol)
        }
        else {
            u_i <- check.KKT(myans)
            newindex <- cpp_findLocalMax(-u_i, main.tol)
        }
        index2 <- sort(unique(c(index, newindex)))
        if (length(index2) - length(index) == 0) {
            k <- length(myans$z)
            myf <- stepfun(z[myans$index], c(0, myans$lam))
            LAMBDA <- myf(z[1:(k - 1)])
            object <- list()
            object$hazard <- myf
            survivalF <- function(x) {
                sapply(x, EstimateSurv, lambda = LAMBDA, z = myans$z)
            }
            class(survivalF) <- "SurvivalFunction"
            object$survival <- survivalF
            object$index <- myans$index
            object$a <- a
            object$b <- b
            object$t <- t
            object$conv <- myans$conv
            object$call <- match.call()
            class(object) <- "iHaz"
            return(object)
        }
        else {
            index <- index2
        }
    }
    k <- length(myans$z)
    myf <- stepfun(z[myans$index], c(0, myans$lam))
    LAMBDA <- myf(z[1:(k - 1)])
    object <- list()
    object$hazard <- myf
    survivalF <- function(x) {
        sapply(x, EstimateSurv, lambda = LAMBDA, z = myans$z)
    }
    class(survivalF) <- "SurvivalFunction"
    object$survival <- survivalF
    object$index <- myans$index
    object$a <- a
    object$b <- b
    object$t <- t
    object$call <- match.call()
    object$conv <- FALSE
    class(object) <- "iHaz"
    return(object)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
