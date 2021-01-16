\name{MonoMat}
\alias{MonoMat}
\title{Matrix of derivative coefficients for local polynomial estimates}
\description{
  This computes a matrix of coefficients of 
the first derivatives of monotonic local linear sharpening problem.
}
\usage{
MonoMat(xgrid, x,  h, d)
}
\arguments{
  \item{xgrid}{numeric vector of locations where monotonicity constraint is to be enforced}
  \item{x}{numeric explanatory vector}
  \item{h}{numeric bandwidth}
  \item{d}{local polynomial degree, can be either 0 or 1}
}
\value{
 a list containing the A matrix and the number of rows in A.  
}
\author{W.J. Braun}
\keyword{models}
