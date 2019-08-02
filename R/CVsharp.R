CVsharp <-
function (x, y, deg, nsteps) 
{
    n <- length(x)
    if (n < 8) 
        stop("Error: not enough observations.  This procedure requires at least 8.")
    y <- y[order(x)]
    x <- sort(x)
    y.lm <- lm(y ~ x + I(x^2))
    gap <- max(diff(sort(x))[-(n-1)]+diff(sort(x))[-1]) # largest distance between pairs of points containing exactly one other point
    hmed <- (summary(y.lm)$sigma^2/(2 * sqrt(pi))/length(x)/(2 * 
        coef(y.lm)[3])^2)^(1/5)
    hmax <- min(hmed * 5, diff(range(x)))
    hmin <- max(hmed/12, gap/3)
    bw <- exp(seq(log(hmin), log(hmax), len = 25))
    CVscore <- numeric(length(bw))
    for (j in 1:length(bw)) {
        resid1 <- NULL
        for (i in 3:(n - 2)) {
            samp <- i
            x1 <- x[-samp]
            y1 <- y[-samp]
            if (nsteps > 0) {
                y1sharp <- sharpiteration(x1, y1, deg, bw[j], 
                  nsteps)[[nsteps]]
            }
            else {
                y1sharp <- y1
            }
            y1.lp <- locpoly(x1, y1sharp, bandwidth = bw[j], 
                degree = deg)
            resid1 <- c(resid1, y[samp] - approx(y1.lp$x, y1.lp$y, 
                x[samp])$y)
        }
        CVscore[j] <- sum(resid1^2)
    }
    CVh <- bw[which.min(CVscore)]
    list(bw = bw, CVscore = CVscore, CVh = CVh)
}
