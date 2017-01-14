# Description: 
#   Plot KL(p,q) and KL(q,p) for 2d Gaussian.
#
# Details:
#   Chapter 21 of Kevin Murphy's book `Machine Learning: a 
#   probabilistic approach` provides a detailed explanation.
#   It doesn't provide much in the way of code though.  
#   This Gist is a brief demo of the KL(p,q) and KL(q,p) 
#   for 2d Gaussian, as described on pages 734.
#
# License: GNU v3
#
# Date: October 18, 2016
#
# Authors:
#    Rcode was modified by Jiali Lin.
#    Depts. of Statistics, Virginia Tech,
#    Hutcheson Hall, 403K, Blacksburg, VA 24061 
#
# References:
#    https://pmtk3.googlecode.com.

rm(list = ls())

library(matlab)
library(mvtnorm)

mu <- c(0, 0)

Sigma <- matrix(c(1, 0.97, 0.97, 1), ncol = 2);
SigmaKLa <- diag(2)/25;
SigmaKLb <- diag(2);


x1 <- seq(-1, 1, 0.01);
x2 <- x1;

n1 <- length(x1);
n2 <- length(x2);

f <- matrix(0, nrow = n1, ncol = n2); 
klqp <- f;
klpq <- f;

for (i in 1:n1){
  f[i, ] <- dmvnorm(cbind(repmat(x1[i], n2, 1), x2), mu, Sigma);
  klqp[i, ] <- dmvnorm(cbind(repmat(x1[i], n2, 1), x2), mu, SigmaKLa);
  klpq[i, ] <- dmvnorm(cbind(repmat(x1[i], n2, 1), x2), mu, SigmaKLb);
}

# figure(a)
contour(x1, x2, f, col = "blue")
contour(x1, x2, klqp, add = TRUE, col = "red") 

# figure(b)
contour(x1, x2, f, col = "blue")
contour(x1, x2, klpq, add = TRUE, col = "red") 