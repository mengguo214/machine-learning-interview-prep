# Description: 
#   This script shows SMC method using SIR algorithm in linear Gaussian state-space model
# 
# Details:
#   x.1 ~ N(0, 1) 
#   x.n = phi x.{n-1} + v.n, v.n ~ N(0,s.v^2) 
#   y.n = x.n + w.n, w.n~ N(0,s.w^2) 
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
#    <http://www.apps.stat.vt.edu/zhu/teaching/2016/6474/6474_2016.htm>

rm(list = ls())

set.seed(100)

# Initial value
T <- 100
S.v <- 1
S.w <- 1 
phi <- 0.95
x <- rep(NA, T)  # latent sequence
y <- rep(NA, T)  # observed sequence
x[1] <- rnorm(1, mean = 0, sd = 1)
for (i in 2:T){
  x[i] <- phi*x[i-1] + rnorm(1, mean = 0, sd = S.v)
}
y <- x + rnorm(T, mean = 0, sd = S.w)

# Initialize storage
N <- 1000 
Xs <- matrix(NA, N, T) 
logweight <- matrix(NA, N, T)
weight.normalized <- matrix(NA, N, T)
mu.Xs <- rep(NA, T)
var.Xs <- rep(NA, T)


for (t in 1:T){
  
  if (t == 1){
    Xs[, t] <- rnorm(N, mean =  0, sd = 1)
  }else{
    # propagate particles according to prior 
    Xs[, t] <- rnorm(N, mean = phi * Xs[, t-1], sd = S.v)    
  }
  logweight[, t] <- dnorm(y[t], mean = Xs[, t], sd = S.w, log = T)    
  
  weight.normalized[, t] <- exp(logweight[, t]) / sum(exp(logweight[, t]))
  
  # Resample
  Xs.resample <- sample(Xs[, t], size = N, replace = T, prob = weight.normalized[, t])
  Xs[, t] <- Xs.resample
  mu.Xs[t] <- mean(Xs.resample)
  var.Xs[t] <- var(Xs.resample)  
}

# win.graph()
plot(1:T, x, type = 'b', pch = 3, col = 'grey', main = 'SIR Solution', ylim = c(-4, 5))
points(1:T, y, pch = 16, col = 'blue')
lines(1:T, y,  col = 'blue')
lines(mu.Xs - 2 * sqrt(var.Xs), lty = 2, col = 'black')
lines(mu.Xs + 2 * sqrt(var.Xs), lty = 2, col = 'black')
lines(mu.Xs, col = "red", lwd = 2)
points(mu.Xs, col = "red", lwd = 2, pch = 16)
legend(10, -2.4, c('latent', 'observed', 'post.mean', 'post.CI'),
       pch = c(3, 16, 16, NA), col = c('grey', 'blue', 'red', 'black'), lty = c(1, 1, 1, 2))