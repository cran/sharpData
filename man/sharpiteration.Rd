\name{sharpiteration}
\alias{sharpiteration}
\title{Iterated Data Sharpening for Local Polynomial Regression}
\description{
Calculation of sharpened responses for bias reduction in function
and estimation, assuming a gaussian kernel is used
in bivariate scatterplot smoothing.
}
\usage{
sharpiteration(x, y, deg, h, nsteps, na.rm, ...)
}
\arguments{
  \item{x}{a numeric vector containing the predictor variable values.}
  \item{y}{a numeric vector containing the response variable values.}
  \item{deg}{a numeric vector containing the local polynomial degree used.}
  \item{h}{a numeric vector containing the (scalar) bandwidth.}
  \item{nsteps}{a numeric vector containing the number of iteration steps.}
  \item{na.rm}{a logical value indicating whether to remove missing values from fitted vectors}
  \item{...}{additional arguments to locpoly}
}
\value{
  a list with elements containing the sharpened (i.e. perturbed) response values, 
ready for input into a local polynomial regression estimator.  The ith list
element corresponds to i steps of data sharpening.
  }
\author{W.J. Braun }
\seealso{ \code{locpoly} }
\examples{
speed <- MPG[, 1]
mpg <- MPG[, 2]
h <- dpill(speed, mpg)
mpgSharp <- sharpiteration(speed, mpg, 1, h, 2)
mpg.lS <- locpoly(speed, mpgSharp[[2]], bandwidth=h, degree=1)
mpg.lX <- locpoly(speed, mpg, bandwidth=h, degree=1)
plot(mpg ~ speed)
lines(mpg.lX)  # unsharpened function estimation
lines(mpg.lS, col=2, lty=2)  # sharpened function estimation (2 steps)
}
\keyword{models}
