sharpiteration <- function(x, y, deg, h, nsteps, na.rm=TRUE, ...) {
    y.sharp <- y
    y.sharp.lst <- vector(nsteps, mode="list")
    sharpstep <- function(x, y, y.sharp, deg, h, na.rm, ...) {
        y.lp <- locpoly(x, y.sharp, degree = deg, bandwidth = h,  gridsize = 201L, ...)
        if (na.rm) {
            cc <- complete.cases(y.lp)
            y.lp$x <- y.lp$x[cc]
            y.lp$y <- y.lp$y[cc]
        }
        g.hat <- approx(y.lp$x, y.lp$y, xout = x)
        y.sharp <- y + y.sharp - g.hat$y
        return(y.sharp)
    }
    for (i in 1:nsteps) {
        y.sharp <- sharpstep(x, y, y.sharp, deg, h, na.rm, ...)
        y.sharp.lst[[i]] <- y.sharp
    }
return(y.sharp.lst)
}

