LLsharpen <- function(x, y, h) {
    ysharp <- y
    for (j in 1:length(x)) {
        y.lm <- lm(y ~ I(x-x[j]), weights=dnorm(x-x[j], sd = h))
        ysharp[j] <- 2*y[j] - coef(y.lm)[1]
    }
    ysharp
}

