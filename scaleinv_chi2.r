dscinvchisq <- function(x,v,s2) {
    v2 <- v/2
    (x*v2)^v2 / gamma(v2) * exp(-v2*x/s2) / s2^(v2+1)
}

x <- seq(0.01,4,length=200)
n <- 5
s2 <- 2
chisq <- dchisq(n*x/s2,n+2)
scinvchisq <- dscinvchisq(x,n,s2)
scinvchisq <- scinvchisq/max(scinvchisq)*max(chisq)

plot(chisq ~ x)
lines(x, scinvchisq, col=2)
