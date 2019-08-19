"NWmono" <-
function(x, y, h, xgrid, numgrid=401, kernel="biweight", call.plot=TRUE, ...){ 
n <- length(y)
if (missing(xgrid)) {xgrid <- seq(range(x)[1], range(x)[2], length=numgrid)}
Dmat <- 2*diag(rep(1, n))  # matrix of quadratic form for objective
Amat <- NWmat(xgrid, x, h, kernel)  # matrix for linear constraints
dvec <- 2*y                # vector for linear form for objective
ysharp <- solve.QP(Dmat, dvec, Amat)$solution
ygrid <- NW.fit(x, ysharp, xgrid, h)
if (call.plot) {plot(x,y,pch=16,...)
points(x, ysharp, col="red")
lines(xgrid,ygrid,col="blue")}
list(x=x,y=y,ysharp=ysharp,h=h,xgrid=xgrid,ygrid=ygrid)
}
