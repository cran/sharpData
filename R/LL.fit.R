"LL.fit" <-
function (x, y, xgrid, h, kernel="biweight")
{
n <- length(x)
xobs <- x
yobs <- y
m <- length(xgrid)
ygrid <- xgrid
work <- x
choice <- seq(1,3)[c("biweight", "normal", "epanechnikov")==kernel]
z <- .Fortran("ll", as.double(xobs), as.double(yobs), as.double(xgrid),
as.double(ygrid), as.integer(n),
as.integer(m), as.integer(choice), as.double(h), as.double(work),
package="sharpData")
names(z) <- c("x", "y", "xgrid", "ygrid", "n", "m", "choice", "h", "work")
z$ygrid
}



