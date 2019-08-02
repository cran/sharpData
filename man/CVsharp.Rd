\name{CVsharp}
\alias{CVsharp}
\title{Cross-Validation Bandwidth Selector for Local Polynomial Regression}
\description{
Cross-validation bandwidth selector for iterated sharpened responses for bias reduction in function
estimation.
}
\usage{
CVsharp(x, y, deg, nsteps)
}
\arguments{
  \item{x}{a numeric vector containing the predictor variable values.}
  \item{y}{a numeric vector containing the response variable values.}
  \item{deg}{a numeric vector containing the local polynomial degree used.}
  \item{nsteps}{a numeric vector containing the number of iteration steps.}
}
\value{
  a list containing 3 elements: the candidate bandwidths; the corresponding
CV scores; the selected optimal bandwidth.
  }
\details{If nsteps is specified to be 0, then the CV bandwidth for conventional
local polynomial regression is provided.
}
\author{W.J. Braun }
\seealso{ \code{locpoly} }
\examples{
speed <- MPG[, 1]
mpg <- MPG[, 2]
h <- CVsharp(speed, mpg, 0, 0)$CVh # conventional local constant regression bandwidth
mpg.l0 <- locpoly(speed, mpg, bandwidth=h, degree=0)
h <- CVsharp(speed, mpg, 0, 1)$CVh # 1-sharpened local constant regression bandwidth
mpgSharp <- sharpiteration(speed, mpg, 0, h, 1)
mpg.l1 <- locpoly(speed, mpgSharp[[1]], bandwidth=h, degree=0)
h <- CVsharp(speed, mpg, 0, 5)$CVh # 5-sharpened local constant regression bandwidth
mpgSharp <- sharpiteration(speed, mpg, 0, h, 5)
mpg.l5 <- locpoly(speed, mpgSharp[[5]], bandwidth=h, degree=0)
plot(mpg ~ speed)
lines(mpg.l0)  # unsharpened function estimation
lines(mpg.l1, col=2, lty=2)  # sharpened function estimation (1 steps)
lines(mpg.l5, col=4, lty=3)  # sharpened function estimation (5 steps)
}
\keyword{models}
