\name{LL.fit}
\alias{LL.fit}
\title{Local linear Regression Estimation}
\description{
  This computes the local linear regression estimator
at a few gridpoints.  For a large number of gridpoints, it is not efficient.
}
\usage{
LL.fit(x, y, xgrid, h, kernel)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{numeric explanatory vector}
  \item{y}{numeric response vector}
  \item{xgrid}{numeric vector}
  \item{h}{numeric bandwidth}
  \item{kernel}{character constant}
}
\value{
 a vector of predictor values  
}
\references{Wand and Jones, 1996, Kernel Smoothing, Chapman and Hall.}
\author{W.J. Braun}
\seealso{\code{\link{locpoly}}}

\examples{
x <- sort(runif(100))
y <- x + sin(x) + rnorm(100)*.1
LL.fit(x,y,c(.3, .5, .7),h=.5)
}
\keyword{models}
