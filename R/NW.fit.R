"NW.fit" <-
function (x,y,xgrid,h,kernel="biweight") 
{
K <- function (x, h, kernel="biweight") {
    if (kernel=="biweight") {
        xh <- x/h
        (abs(xh) < 1)*(1-xh^2)^2*15/(16*h)  # biweight
    } else {
        dnorm(x, sd=h)
    }
}
subtract <- function(x,y){x-y}
diffs <- outer(xgrid, x, subtract)
denominator <- K(diffs, h, kernel)%*%rep(1,length(y))
numerator <- K(diffs, h, kernel)%*%y
numerator/denominator
}
