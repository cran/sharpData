MonoMat <- function (xgrid, x, h, d) 
{
    n <- length(x)
    xobs <- x
    m <- length(xgrid)
    amat <- rep(0, n * m)
    loclin <- dloclin <- numeric(n)
    z <- .Fortran(C_afun, as.double(xobs), as.double(xgrid), as.integer(n), 
        as.integer(m), as.double(amat), as.double(loclin), as.double(dloclin), 
        as.double(h), as.integer(d))
    names(z) <- c("x", "xgrid", "n", "m", "Amat", "loclin", "dloclin", 
        "h", "d")
    z$Amat <- matrix(z$Amat, nrow = n)
    if (d == 1) z$Amat <- z$Amat[, 1:(m-1)]
    z$Amat
}
