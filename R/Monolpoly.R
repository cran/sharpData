"Monolpoly" <- function(x, y, h, d = 1, xgrid, numgrid=401, ...){ 
    n <- length(y)
    if (missing(xgrid)) {
        xgrid <- seq(range(x)[1L], range(x)[2L], length=numgrid)
    }
    Dmat <- 2*diag(rep(1, n)) # matrix of quadratic form for objective if (d == 1) {
    Amat <- MonoMat(xgrid, x, h, d)  # matrix for NW constraints
    dvec <- 2*y                # vector for linear form for objective
    ysharp <- solve.QP(Dmat, dvec, Amat)$solution
    sharp.out <- locpoly(x, ysharp, bandwidth=h, degree=d, ...)
    list(x=sharp.out$x, y=sharp.out$y, ysharp=ysharp)
}
