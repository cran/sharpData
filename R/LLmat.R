LLmat <- function (xgrid, x, h, kernel = "biweight") 
{
    n <- length(x)
    xobs <- x
    m <- length(xgrid)
    amat <- rep(0, n * m)
    loclin <- dloclin <- work <- numeric(n)
    choice <- seq(1, 3)[c("biweight", "normal", "epanechnikov") == 
        kernel]
    z <- .Fortran("dllc", as.double(xobs), as.double(xgrid), as.integer(n), 
        as.integer(m), as.double(amat), as.double(loclin), as.double(dloclin), 
        as.integer(choice), as.double(h), as.double(work), PACKAGE = "sharpData")
    names(z) <- c("x", "xgrid", "n", "m", "Amat", "loclin", "dloclin", 
        "choice", "h", "work")
    matrix(z$Amat, nrow = n)
}
