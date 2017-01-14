# Description: 
#   This script shows rejection sampling from target density.
#
# Details:
#   Target density: sqrt(2/pi)*exp(-theta.list^2/2) 
#   Proposal distribution: g.theta = exp(-theta.list)
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
# References:
#    http://www.apps.stat.vt.edu/zhu/teaching/2016/6474/6474_2016.htm
# 
rm(list = ls())

theta.list <- seq(0, 10, length.out = 100)
f.theta <- sqrt(2 / pi) * exp(-theta.list^2 / 2)
g.theta <- exp(-theta.list)
M <- 1.4

plot(theta.list, f.theta, type = "l", lwd = 2)
lines(theta.list, g.theta*M, col = "blue", lwd = 3)

sample <- rep(NA, 1000)
tag <- 0
while (tag <= 1000){
  theta.i <- rexp(1, rate = 1)
  f.thetai <- sqrt(2 / pi) * exp(-theta.i^2 / 2)
  M.gthetai <- M*dexp(theta.i, rate = 1)
  prob.i <- f.thetai/M.gthetai
  u <- runif(1)
  if (u < prob.i){
    tag <- tag+1
    sample[tag] <- theta.i 
  }
}

hist(sample, breaks = 20, probability = TRUE, add = TRUE, col = "orange")
