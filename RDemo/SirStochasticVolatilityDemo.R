# Description: 
#   This script shows SMC method using SIR algorithm in Stochastic volatility model
#
# Details:
#   x.1 ~ N(0, \frac{\sigma^2}{1-\sigma^2})
#   x.n = phi x.{n-1} + v.n, v.n ~ N(0,\sigma^2), |sigma^2|<1
#   y.n = exp{gamma + x.n}w.n, w.n~ N(0,1)
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

# Initial Value
Time <- 100
sigma <- 1
gamma <- 1 
phi <- 0.9
x <- rep(NA, Time)  # latent sequence
y <- rep(NA, Time)  # observed sequence
x[1] <- rnorm(1, mean = 0, sd = sigma / sqrt(1 - phi^2))
for (i in 2:Time){
  x[i] <- phi * x[i-1] + rnorm(1, mean = 0, sd = sigma)
}
y <- exp(gamma + x) * rnorm(Time, mean = 0, sd = 1)

# Plot the latent sequence and the observed sequence
# win.graph()
plot(1:Time, x, type = 'b', pch = 3, col = 'green', xlab = "Time", ylab = NULL)
points(1:Time, y, pch = 16, col = 'blue')
lines(1:Time, y,  col = 'blue')
legend(-1, -3, c('latent', 'observed'), col = c('green','blue'), pch = c(3, 16), lty = c(1, 1))


#  SMC method using SIR algorithm.
# -------------
N <- 1000 
Xs <- matrix(NA, N, Time)  
logweight <- matrix(NA, N, Time) 
weight.normalized <- matrix(NA, N, Time)
mu.Xs <- rep(NA, Time)
var.Xs <- rep(NA, Time)

# sample from p(x_n|y_{1:n})
for (t in 1:Time){
  
  if (t == 1){
    # Xs[,t] = rnorm(N, mean= 0, sd = sigma/sqrt(1+phi^2)) 
    Xs[, t] <- rnorm(N, mean = 0, sd = sigma / sqrt(1 - phi^2)) #**** zhunote: 1-phi^2
  }else{
    # propagate particles according to prior 
    Xs[, t] <- rnorm(N, mean = phi * Xs[, t-1], sd = sigma)    
  }
  logweight[, t] <- dnorm(y[t], mean = 0, sd = exp(gamma + Xs[, t]), log = T)    
  
  weight.normalized[,t] <- exp(logweight[, t])/sum(exp(logweight[, t]))
  
  # Resample
  Xs.resample <- sample(Xs[, t], size = N, replace = T, prob = weight.normalized[, t])
  Xs[, t] <- Xs.resample
  mu.Xs[t] <- mean(Xs.resample)
  var.Xs[t] <- var(Xs.resample) 
}

# win.graph()
plot(1:Time, x, type = 'b', pch = 3, col = 'grey', main = 'SIR Solution', ylim = c(-10, 20), xlab = "Time")
points(1:Time, y, pch = 16, col = 'blue')
lines(1:Time, y,  col = 'blue')
lines(mu.Xs - 2 * sqrt(var.Xs), lty = 2, col = 'black') 
lines(mu.Xs + 2 * sqrt(var.Xs), lty = 2, col = 'black')
lines(mu.Xs, col = "red",lwd = 2)
points(mu.Xs, col = "red",lwd = 2, pch = 16)

legend(0, -5, c('latent', 'observed', 'post.mean', 'post.CI'),
       pch = c(3, 16, 16, NA), col = c('grey', 'blue', 'red', 'black'), lty = c(1, 1, 1, 2))