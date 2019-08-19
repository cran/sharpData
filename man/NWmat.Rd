\name{NWmat}
\alias{NWmat}
\title{Smoother Matrix for Nadaraya-Watson Local Constant Regression}
\description{
  This function computes the A matrix for use in quadprog.  This
matrix corresponds to the derivative matrix for constrained
sharpening in the monotonic regression problem.
}
\usage{
NWmat(xgrid, x, h, kernel = "biweight")
}
\arguments{
  \item{xgrid}{a vector of gridpoints to evaluate estimate at}
  \item{x}{a vector containing the predictor data}
  \item{h}{a scalar bandwidth}
  \item{kernel}{a character string naming the kernel function}
}
\value{
  a matrix
}
\author{W.J. Braun}
\keyword{models}
