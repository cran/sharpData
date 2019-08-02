\name{LLsharpen}
\alias{LLsharpen}
\title{Data Sharpening for Local Linear Regression}
\description{
Calculation of sharpened responses for bias reduction in function
and first derivative estimation, assuming a gaussian kernel is used
in bivariate scatterplot smoothing.
}
\usage{
LLsharpen(x, y, h)
}
\arguments{
  \item{x}{a numeric vector containing the predictor variable values.}
  \item{y}{a numeric vector containing the response variable values.}
  \item{h}{a numeric vector containing the (scalar) bandwidth.}
}
\value{
  a vector containing the sharpened (i.e. perturbed) response values, 
ready for input into a local linear regression estimator.
  }
\author{W.J. Braun }
\seealso{ \code{locpoly} }
\examples{
speed <- MPG[, 1]
mpg <- MPG[, 2]
h <- dpill(speed, mpg)*2
mpgSharp <- LLsharpen(speed, mpg, h)
mpg.lS <- locpoly(speed, mpgSharp, bandwidth=h, drv=1, degree=1)
mpg.lX <- locpoly(speed, mpg, bandwidth=h, drv=1, degree=1)
plot(mpg.lX, type="l")  # unsharpened derivative estimation
lines(mpg.lS, col=2, lty=2)  # sharpened derivative estimation
}
\keyword{models}
