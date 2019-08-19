\name{LLmono}
\alias{LLmono}
\title{Monotonized Local Linear Regression}
\description{
  Locally linear regression is applied to bivariate
data.  The response is `sharpened' or perturbed in a way to render a
monotonically increasing curve estimate.
}
\usage{
LLmono(x, y, h, xgrid, numgrid = 401, kernel="biweight", call.plot = 
TRUE, ...)
}
\arguments{
  \item{x}{a vector of explanatory variable observations}
  \item{y}{binary vector of responses}
  \item{h}{bandwidth}
  \item{xgrid}{gridpoints on x-axis where estimates are taken}
  \item{numgrid}{number of equally-spaced gridpoints (if xgrid not 
specified)}
  \item{kernel}{character constant}
  \item{call.plot}{if TRUE (default), the data, sharpened data and 
estimated curve
are plotted.}
  \item{...}{other arguments for plot}
}
\details{
  Data are perturbed the smallest possible L2 distance subject to the 
constraint
that the local linear estimate is monotonically increasing.
}
\value{
  \item{x}{original explanatory variable}
  \item{y}{original responses}
  \item{ysharp}{sharpened responses}
  \item{h}{bandwidth}
  \item{xgrid}{gridpoints}
  \item{ygrid}{sharpened curve estimate}
}
\references{Braun, W.J. and Hall, P., Data Sharpening for Nonparametric 
Estimation
Subject to Constraints, Journal of Computational and Graphical Statistics, 
2001}
\author{W.J.Braun}
\examples{
gridpts <- seq(1, 10, length=101)
x <- seq(1, 10, length=51)
p <- exp(-1 + .2*x)/(1 + exp(-1 + .2*x))
y <- rbinom(51, 1, p)
LLmono(x, y, h=0.6, xgrid=gridpts)
lines(x,p) # true mean response
}
\keyword{models}
