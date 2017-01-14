# Description: 
#   This script shows sampler for univariate random variate using its CDF
#
# Details:
#   Goal: sample from   f = ((x-1)^2) * exp(-(x^3/3-2*x^2/2+x))
#   f(x)     = density function, 
#   F(x)     = cdf of f(x), and 
#   F.inv(y) = inverse cdf of f(x).
# 
# License: GNU v3
#
# Date: October 19, 2016
#
# Authors:
#    Rcode was written by Jiali Lin.
#    Depts. of Statistics, Virginia Tech,
#    Hutcheson Hall, 403K, Blacksburg, VA 24061 
#
# Source:
#    http://stackoverflow.com/questions/20508400/generating-random-sample-from-the-quantiles-of-unknown-density-in-r

rm(list = ls())

f <- function(x) {((x - 1)^2) * exp(-(x^3 / 3 - 2 * x^2 / 2 + x))}
F <- function(x) {integrate(f, 0, x)$value}
F <- Vectorize(F)

F.inv <- function(y){uniroot(function(x){F(x) - y}, interval = c(0, 10))$root}
F.inv <- Vectorize(F.inv)

x <- seq(0, 5, length.out = 1000)
y <- seq(0, 1, length.out = 1000)

par(mfrow=c(1, 3))
plot(x, f(x), type = "l", main = "f(x)")
plot(x, F(x), type = "l", main = "CDF of f(x)")
plot(y, F.inv(y), type = "l", main = "Inverse CDF of f(x)")

# generate a random sample 
X <- runif(1000, 0, 1)  
Z <- F.inv(X)

par(mfrow=c(1,2))
plot(x,f(x),type="l",main="Density function")
hist(Z, breaks=20, xlim=c(0,5))
