\name{MPG}
\alias{MPG}
\title{Mileage Data}
\usage{data(MPG)}
\description{
The \code{MPG} data frame has 15 rows and 10 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{speed}{a numeric vector of cruising speeds in miles per hour}
    \item{corsica88}{miles per gallon for a 1988 Corsica}
    \item{legacy93}{miles per gallon for a 1993 Legacy}
    \item{olds94}{miles per gallon for a 1994 Oldsmobile}
    \item{cutlass94}{miles per gallon for a 1994 Oldsmobile Cutlass}
    \item{chevpickup94}{miles per gallon for a 1994 Chevrolet Pickup}
    \item{cherokee94}{miles per gallon for a 1994 Jeep Cherokee}
    \item{villager94}{miles per gallon for a 1994 Villager}
    \item{prizm95}{miles per gallon for a 1995 Prizm}
    \item{celica97}{miles per gallon for a 1997 Toyota Celica}
  }
}
\source{
B.H.  West, R.N.  McGill, J.W.  Hodgson, S.S.  Sluder, D.E.  Smith,
Development and Verification of Light-Duty Modal Emissions and Fuel
Consumption Values for Traffic Models, Washington, DC, April 1997, and
additional project data, April 1998.
}
\examples{
data(MPG)
plot(celica97 ~ speed, data = MPG)
}
\keyword{datasets}
