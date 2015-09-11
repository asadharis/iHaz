\name{plot.iHaz}
\alias{plot.iHaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots of the estimated hazard and survival functions
}
\description{
Plot function for objects of class "iHaz". This function plots the estimated hazard function and the estimated survival function. This is a generic \code{S3} plot function and accepts additional arguments such as line type, color etc.
}
\usage{
\method{plot}{iHaz}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 An object of class "iHaz", the output of our \code{iHaz} function.
}
  \item{\dots}{
 Additonal arguments to be passed to the plot function. All additional arguments can be passed to the plot functions except \code{main, xlab} and \code{ylab}.
}
}
\author{
Asad Haris, Gary Chan
}

\seealso{
\code{\link{iHaz}}
}
\examples{
library(iHaz)
#Generate some data
#Here event time is distributed Exponential(1)
#Simulation study from Pan et al. (1998)
n<- 500
x<- rexp(n)
t<- runif(n, min = 0, max = 1.5)
xnew<- x[(x>=t)]
tnew<- t[x>=t]
t<- tnew
a<- xnew
b<- xnew
a[xnew<= tnew+0.5]<- tnew[xnew<= tnew+0.5]
b[xnew<= tnew+0.5]<- tnew[xnew<= tnew+0.5] +0.5
a[xnew > tnew+0.5]<- tnew[xnew > tnew+0.5] +0.5
b[xnew > tnew+0.5]<- Inf
dat<- list("a" = a, "b" = b, "t" = t)

#Fit an 'iHaz' object
fit<- iHaz(dat, ini.index = 1:3 ,verbose =TRUE)

#Veiw/plot extimated hazard function
fit$hazard
plot(fit$hazard, main = "Hazard Function")

#View the survival function
fit$survival
#Estimated survival probabilitcies at some time points
fit$survival(c(0.5,0.8,1,1.5))

#Code NOT RUN
#Plot iHaz object
#plot(fit, col = "red", type = "o", lwd = 1, pch = 16, cex = 0.5)
}


