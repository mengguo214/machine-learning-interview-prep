# Description: 
#   This script shows an example of stick breaking procss:
#   Assume F ~ DP(alpha, G0 with alpha = 5, G0 ~ N (0, 1)$. 
#   Approximating F by F.hat = sum (w.h * delta(theta.h)), where H = 100
# 
# Details:
#   Return a vector of weights drawn from a stick-breaking process with dispersion `α`.
#   
#   Recall that the kth weight is
#   beta.k = (1 - beta.1) * (1 - beta.2) * ... * (1 - beta.{k-1}) * beta.k
#     where each `beta.i` is drawn from a Beta distribution
#   beta.i ~ Beta(1, α)
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
#    http://www.apps.stat.vt.edu/zhu/teaching/2016/6474/6474_2016.htm

rm(list = ls())

# Version 1
# ------------
# generates from stick-breaking construction
alpha <- c(1, 5, 10);
nn <- 59;

par(mfrow = c(length((alpha)), 1)) 

for(ii in 1:length(alpha)){
  beta <- rbeta(nn, 1, alpha[ii]);
  neg <- cumprod(1 - beta);
  pi <- beta * c(1, neg[1:length(neg) - 1]);
  hist(pi, main = 'stick-breaking weights pi', ylab = 'alpha', col = "black")
}

# Version 2
# ------------
alpha <- 5
G_0 <- function(n) rnorm(n, 0, 1)
H <- 10
V <- rbeta(H, 1, alpha)
w <- numeric(H)
w[1] <- V[1]
w[2:H] <- sapply(2:H, function(i) V[i] * prod(1 - V[1:(i-1)]))
theta_star <- G_0(H)
theta <- sample(theta_star, prob = w, replace = TRUE)
# Plot $F.hat$. 
par(mfrow = c(1, 1)) 
plot(theta_star, w, type="h") 
