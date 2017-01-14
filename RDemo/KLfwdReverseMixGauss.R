# Description: 
#   Plot KL(p,q) and KL(q,p) for a mix of 2d Gaussian.
#
# Details:
#   Visualize difference between KL(p,q) and KL(q,p) where p 
#   is a mix of two 2d Gaussians, and q is a single 2d Gaussian.
#   This is Figure bishop-5-4.
#
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

mu <- matrix(c(-1, -1, 1, 1), nrow = 2,byrow = T);

Sigma <- array(0,c(2, 2, 2));
Sigma[, , 1] <- matrix(c(1/2, 1/4, 1/4, 1), nrow = 2);
Sigma[, , 2] <- matrix(c(1/2, -1/4, -1/4, 1), nrow = 2);
SigmaKL <- matrix(c(3, 2, 2, 3), nrow = 2);

x1 <- seq(-4, 4, 0.1);
x2 <- x1;

n1 <- length(x1);
n2 <- length(x2);

f1 <- matrix(0, nrow = n1, ncol = n2); 
f2 <- matrix(0, nrow = n1, ncol = n2); 
klf <- matrix(0, nrow = n1, ncol = n2); 
kll <- klf;
klr <- klf;

for (i in 1:n1){
  f1[i,] <- dmvnorm(cbind(repmat(x1[i], n2, 1), x2), mu[1, ], Sigma[, , 1]);
  f2[i,] <- dmvnorm(cbind(repmat(x1[i], n2, 1), x2), mu[2, ], Sigma[, , 2]);
  klf[i,] <- dmvnorm(cbind(repmat(x1[i], n2, 1), x2), matrix(c(0, 0), ncol = 1), SigmaKL);
  kll[i,] <- dmvnorm(cbind(repmat(x1[i], n2, 1), x2), mu[1, ], Sigma[, , 1] * 0.6);
  klr[i,] <- dmvnorm(cbind(repmat(x1[i], n2, 1), x2), mu[2, ], Sigma[, , 2]*  0.6);
}

f <- f1 + f2;

# figure(a)
contour(x1, x2, f, col = "blue")
contour(x1, x2, klf, add = TRUE, col = "red") 

# figure(b)
contour(x1, x2, f, col = "blue")
contour(x1, x2, kll, add = TRUE, col = "red") 

# figure(b)
contour(x1, x2, f, col = "blue")
contour(x1, x2, klr, add = TRUE, col = "red") 