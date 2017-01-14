# Description: 
#   This script shows Metropolis Hasting samlinging MCMC for 1-d Gaussian Mixture Model.
# 
# License: GNU v3
#
# Date: October 19, 2016
#
# Authors:
#    Rcode was written by Jiali Lin.
#    Depts. of Statistics, Virginia Tech,
#    Hutcheson Hall, 403K, Blacksburg, VA 24061 

rm(list = ls())

# Target and proposal distribution 
p <- function(x) 0.4*dnorm(x,0,1) + .6*dnorm(x,5,1)
q <- function(mu) rnorm(1,2,1)                     

# initial
N <- 10000              
xs <- rep(0,N)       
xs[1] <- 0.3               

# main algorithm
for (i in 1:N){
  xs_ast <- q()          
  alpha <- p(xs_ast) / p(xs[i])  
  
  # Accept the sample with prob = min(alpha,1)
  if (runif(1) < min(alpha, 1)){
    xs[i+1] <- xs_ast
  } else {
    xs[i+1] <- xs[i]
  }
  
}  

# Density plot, sample histogram, traceplot
par(mfrow = c(1, 3))
x <- seq(-5, 10, .1)
plot(x, p(x), type="l", lwd=2)
#hist(xs, breaks = 20, probability = TRUE, add = TRUE, col = "orange")
hist(xs, breaks = 20, probability = TRUE, col = "orange")
plot(xs, type = "l", xlab ="", ylab = "", xaxt = "n",
     yaxt = "n", bty = "l", col = "blue")
