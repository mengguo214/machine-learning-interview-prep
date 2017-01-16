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
#    https://pmtk3.googlecode.com.

rm(list = ls())

# ------------
# generates from stick-breaking construction
alpha <- c(2, 5)
nn <- 20;

par(mfrow = c(length((alpha)), 2)) 

for(ii in 1:length(alpha)){
  for (trial in 1:2){
    beta <- rbeta(nn, 1, alpha[ii]);
    neg <- cumprod(1 - beta);
    pi <- beta * c(1, neg[1:length(neg) - 1]);
    barplot(pi, main = bquote(~alpha ~ "=" ~.(alpha[ii])), col = "black")
  }
}