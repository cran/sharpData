\name{LLmat}
\alias{LLmat}
\title{Local Linear Regression Matrix}
\description{
  This computes the smoother matrix for the monotonic local linear sharpening problem.
}
\usage{
LLmat(xgrid, x,  h, kernel)
}
\arguments{
  \item{xgrid}{numeric vector}
  \item{x}{numeric explanatory vector}
  \item{h}{numeric bandwidth}
  \item{kernel}{character constant}
}
\value{
 a list containing the A matrix and the number of rows in A.  
}
\author{W.J. Braun}
\keyword{models}
