# Description: 
#   This script shows slice sampling univariate distribution.
# 
# Details:
#   Target Distribution: f(theta) = (1+sin^2(3theta))(1+cos^4(5theta))\exp{-theta^2/2}   
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

nMC <- 2000
theta.sample <- rep(NA, nMC)
u.sample <- matrix(NA,nrow = nMC,ncol = 3)
theta <- 0.5
for (i in 1:nMC){
  # step 1 sample u1, u2, u3 from uniform
  u1 <- runif(1, min = 0, max = 1+(sin(3*theta))^2)
  u2 <- runif(1, min = 0, max = 1+(cos(5*theta))^4)
  u3 <- runif(1, min = 0, max =exp(-theta^2/2))
  u.sample[i,] = c(u1,u2,u3)
  
  # step 2 sample theta from A(u1, u2, u3)
  repeat {
    theta <- runif(1, min = -sqrt(-2 * log(u3)), max = sqrt(-2 * log(u3)))
    if ((sin(3 * theta))^2>= u1 - 1 && (cos(5 * theta))^4 >= u2 - 1) break  
  }
  theta.sample[i] <- theta
}

# win.graph()
burnin <- 1000
p1 <- hist(theta.sample[(burnin+1):nMC], nclass = 75, col = "orange",
          proba = T, xlab = "x", ylab = "", main = "")
theta_grid <- seq(-3, 3, length.out = 100)
f.theta <- (1 + sin(3 * theta_grid)^2) * (1 + cos(5 * theta_grid)^4) * exp(-theta_grid^2 / 2)
f.theta <- f.theta * max(p1$density) / max(f.theta)
lines(theta_grid, f.theta, col = "black", lwd = 2)