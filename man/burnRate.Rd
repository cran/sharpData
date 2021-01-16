\name{burnRate}
\alias{burnRate}
\title{Firebrand Burning Properties}
\usage{data(burnRate)}
\description{
The \code{burnRate} data frame contains laboratory data on the
proportion of remaining fuel in a piece of wood that has burned
for a fixed period of time subjected to a fixed windspeed.  
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{proportionBurned}{a numeric vector}
    \item{densityRatio}{ratio of windspeed, multiplied by density of air,
 to density of firebrand}
    \item{species}{factor listing tree species}
    \item{diameter}{numeric vector of diameter of burned particle in cm}
    \item{windspeed}{windspeed in cm per second}
    \item{testTime}{length of test in seconds}
  }
}
\source{
Albini, F. USDA Forest Service General Technical Report INT-56, 1979.  
}
\keyword{datasets}
