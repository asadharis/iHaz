\name{iHaz-package}
\alias{iHaz-package}

\docType{package}
\title{
Increasing hazard function estimation for truncated and censored data
}
\description{
Increasing hazard function estimation for truncated and interval censored data. This package builds upon the work on Pan et al. (1998) who considered the problem of estimating monotone hazards for left truncated, interval censored data. We incorporate ideas from the support reduction algorithm [Groeneboom et al. (2008)] to efficiently handle larger datasets.
}
\details{

\tabular{ll}{
Package: \tab iHaz\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2015-09-09\cr
License: \tab GPL (>= 2) \cr
}

The package includes the following functions:

\tabular{ll}{
\code{\link{iHaz}}: \tab Estimate the hazard and survival functions \cr
\code{\link{plot.iHaz}}: \tab S3 \code{plot} method for class \code{"iHaz"}. Plots the
estimated hazard and survival function\cr
}
}

\author{
Asad Haris, Gary Chan

Maintainer: Asad Haris <aharis@uw.edu>
}

\references{
Pan, Wei, and R. Chappell. "Estimating survival curves with left-truncated and interval-censored data under monotone hazards." Biometrics (1998): 1053-1060.

Piet Groeneboom, Geurt Jongbloed, and Jon A. Wellner. "The Support Reduction Algorithm for Computing Non-Parametric Function Estimates in Mixture Models." Scandinavian Journal of Statistics 35.3 (2008): 385-399.
}

\keyword{ package }
\seealso{
\code{\link{iHaz}},
\code{\link{plot.iHaz}},
}

