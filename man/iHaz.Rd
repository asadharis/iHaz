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
A string to specify the type of convergence criteria. The possible options are "KKT" which uses the Karush–Kuhn–Tucker conditions, or "derv" which uses a derivative condition as in Lemma 3.1 of Wellner et al. (1997).
}
  \item{verbose}{
Logical indicator to specify whether or not details should be printed to the screen while running the algorithm.
}
}
\details{
The estimated hazard function via the projection algorithm of Pan et al. (1998) is a step function with a potential change point at every time point in our dataset. In practice however, the estimated hazard has only a small number of change points. This naturally leads us to a support reduction type algorithm [Groeneboom et al. (2008)]. We define \code{z = sort(unique(c(a,b,t)))} as all the possible time points. Our estimated hazard function is a step function with potential change points at each element of \code{z}. For then algorithm of Pan et al. (1998) this gives us a parameter vector \code{lambda} of length \code{K = length(z)}. In practice however, the total number of change points for the estimated step function is much smaller than \code{K}.

Our algorithm proceeds as follows:

We begin with \code{ini.index}, the initial indices which are a subset of the vector \code{1:K}. For comupational reasons we must have 1 in the initial index set. We then optimizate the objective function over the class of step functions with change points given by \code{ini.index}.This is done via the projection algorithm of Pan et al. (1998). \code{inter.maxiter} and \code{inter.tol} specify the maximum number of iterations and tolerance for checking stopping conditions for this projection algorithm, respectively.

Once the function is optimized for the reduced support, we check if the some optimization condition is met upto the \code{main.tol} accuracy. If the condition is not met, support points are added and the index set is updated. This continues until convergence or until \code{main.maxiter} iterations.
}
\value{
An object of class "iHaz" whcih is a list with the following components:
\item{hazard}{
The estimated hazard function. This is an \code{R}, vector of class "stepfun".
}
\item{survival}{
The etimated survival function. This function also takes as input, a vector at which we wish to evaluate the estimated survival probability.
}
\item{index}{
The index set of our estimates. This corresponds to the indices of \code{z}; it is the final updated index set of our algorithm.
}
\item{a, b, t}{
The values of \code{a, b} and \code{t} used in estimation of survival and hazard functions.
}
\item{conv}{

}
}
\references{
Pan, Wei, and R. Chappell. "Estimating survival curves with left-truncated and interval-censored data under monotone hazards." Biometrics (1998): 1053-1060.

Piet Groeneboom, Geurt Jongbloed, and Jon A. Wellner. "The Support Reduction Algorithm for Computing Non-Parametric Function Estimates in Mixture Models." Scandinavian Journal of Statistics 35.3 (2008): 385-399.

Jon A. Wellner, and Yihui Zhan. "A hybrid algorithm for computation of the nonparametric maximum likelihood estimator from censored data." Journal of the American Statistical Association 92.439 (1997): 945-959.
}
\author{
Asad Haris, Gary Chan
}

\seealso{
\code{\link{plot.iHaz}}
}
\examples{

}
