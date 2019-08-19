"NWmat" <-
function (xgrid, x, h, kernel="biweight") 
{
n <- length(x)
xobs <- x
m <- length(xgrid)
amat <- rep(0,n*m)
choice <- seq(1,3)[c("biweight", "normal", "epanechnikov")==kernel]
z <- .Fortran("afun", as.double(xobs), as.double(xgrid), as.integer(n),
as.integer(m), as.double(amat), as.integer(choice), as.double(h), PACKAGE="sharpData")
names(z) <- c("x", "xgrid", "n", "m", "Amat", "choice", "h")
matrix(z$Amat, nrow=n)
}


