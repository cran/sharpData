\name{Monolpoly}
\alias{Monolpoly}
\title{Monotonized Local Regression}
\description{
  Local constant and local linear regression are applied to bivariate
data.  The response is `sharpened' or perturbed in a way to render a
monotonically increasing curve estimate. 
}
\usage{
Monolpoly(x, y, h, d=1,  xgrid, numgrid = 401, ...)
}
\arguments{
  \item{x}{a vector of explanatory variable observations}
  \item{y}{binary vector of responses}
  \item{h}{bandwidth}
  \item{d}{degree, can be either 0 or 1}
  \item{xgrid}{gridpoints on x-axis where monotonicity constraint is enforced}
  \item{numgrid}{number of equally-spaced gridpoints (if xgrid not 
specified)}
  \item{...}{other arguments for locpoly}
}
\details{
  Data are perturbed the smallest possible L2 distance subject to the 
constraint
that the local linear estimate is monotonically increasing.
}
\value{
  \item{x}{locations of function estimate evaluations}
  \item{y}{function estimate evaluations (sharpened - monotonized)}
  \item{ysharp}{sharpened responses}
}
\references{Braun, W.J. and Hall, P., Data Sharpening for Nonparametric 
Estimation Subject to Constraints, Journal of Computational and Graphical Statistics, 
2001}
\author{W.J.Braun}
\examples{
gridpts <- seq(1, 10, length=101)
x <- seq(1, 10, length=51)
p <- exp(-1 + .2*x)/(1 + exp(-1 + .2*x))
y <- rbinom(51, 1, p)
plot(x, y)
lines(Monolpoly(x, y, h=0.6, xgrid=gridpts))
##
plot(faithful)
with(faithful, 
lines(Monolpoly(eruptions, waiting, h=0.1, d=1, 
range=c(1.55,5.15))))
}
\keyword{models}
