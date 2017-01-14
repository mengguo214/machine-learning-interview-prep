# Description: 
#   Demo of Gaussian process regression with R.
#
# Details:
#   Chapter 2 of Rasmussen and Williams's book `Gaussian Processes
#   for Machine Learning' provides a detailed explanation of the
#   math for Gaussian process regression.  It doesn't provide
#   much in the way of code though.  This Gist is a brief demo
#   of the basic elements of Gaussian process regression, as
#   described on pages 13 to 16.
#
# License: GNU v3
#
# Date: October 18, 2016
#
# Authors:
#    Originally called GP.R.
#    Rcode was written by James Keirstead.
#    Rcode was modified by Jiali Lin.
#    Depts. of Statistics, Virginia Tech,
#    Hutcheson Hall, 403K, Blacksburg, VA 24061 
#
# References:
#    https://gist.github.com/jkeirstead/2312411
# 
# See Also:
#    http://www.jameskeirstead.ca/blog/gaussian-process-regression-with-r/

# Load in the required libraries for data manipulation
# and multivariate normal distribution

# Set a seed for repeatable plots
rm(list = ls())
set.seed(12345)
require(MASS)

calcSigma <- function(X1, X2, l = 1) {
  
  # Calculates the covariance matrix sigma using a
  # simplified version of the squared exponential function.
  #
  # Parameters:
  #  X1, X2 = vectors
  # 	l = the scale length parameter
  # Returns:
  # 	a covariance matrix
  
  Sigma <- matrix(rep(0, length(X1) * length(X2)), nrow = length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i, j] <- exp(-0.5 * (abs(X1[i] - X2[j]) / l)^2)
    }
  }
  return(Sigma)
}

# 1. Plot some sample functions from the Gaussian process
# as shown in Figure 2.2(a)

# Define the points at which we want to define the functions
x.star <- seq(-5, 5, len = 50)

# Calculate the covariance matrix
sigma <- calcSigma(x.star, x.star)

# Generate a number of functions from the process
n.samples <- 3
values <- matrix(NA, nrow = length(x.star), ncol = n.samples)
for (i in 1:n.samples) {
  # Each column represents a sample from a multivariate normal distribution
  # with zero mean and covariance sigma
  values[, i] <- mvrnorm(1, rep(0, length(x.star)), sigma)
}

plot(x.star, values[, 1], type = 'l', ylim = c(-2, 2))
for (i in 2:n.samples){
  lines(x.star, values[, i])
}


# 2. Now let's assume that we have some known data points;
# this is the case of Figure 2.2(b). In the book, the notation 'f'
# is used for f$y below.  I've done this to make the ggplot code
# easier later on.
f <- data.frame(x = c(-4, -3, -1, 0, 2),
               y = c(-2, 0, 1, 2, -1))

# Calculate the covariance matrices
# using the same x.star values as above
x <- f$x
k.xx <- calcSigma(x, x)
k.xxs <- calcSigma(x, x.star)
k.xsx <- calcSigma(x.star, x)
k.xsxs <- calcSigma(x.star, x.star)

# These matrix calculations correspond to equation (2.19)
# in the book.
f.star.bar <- k.xsx %*% solve(k.xx) %*% f$y
cov.f.star <- k.xsxs - k.xsx %*% solve(k.xx) %*% k.xxs

# This time we'll plot more samples.  We could of course
# simply plot a +/- 2 standard deviation confidence interval
# as in the book but I want to show the samples explicitly here.
n.samples <- 50
values <- matrix(rep(0, length(x.star) * n.samples), ncol = n.samples)
for (i in 1:n.samples) {
  values[, i] <- mvrnorm(1, f.star.bar, cov.f.star)
}

plot(x.star, values[, 1], type = 'l', ylim = c(-3, 3), col = 'gray')
for (i in 2:n.samples){
  lines(x.star, values[, i], col = 'gray')
}
points(f$x, f$y, pch=21, col = 'blue')


# 3. Now assume that each of the observed data points have some
# normally-distributed noise.

# The standard deviation of the noise
sigma.n <- 0.1

# Recalculate the mean and covariance functions
f.bar.star <- k.xsx %*% solve(k.xx + sigma.n^2 * diag(1, ncol(k.xx))) %*% f$y
cov.f.star <- k.xsxs - k.xsx %*% solve(k.xx + sigma.n^2 * diag(1, ncol(k.xx))) %*% k.xxs

# Recalulate the sample functions
values <- matrix(rep(0, length(x.star) * n.samples), ncol = n.samples)
for (i in 1:n.samples) {
  values[, i] <- mvrnorm(1, f.bar.star, cov.f.star)
}

plot(x.star, values[, 1], type = 'l', ylim = c(-3, 3), col = 'gray')
for (i in 2:n.samples){
  lines(x.star, values[, i], col = 'gray')
}
points(f$x, f$y, pch = 21, col = 'blue')
for (i in 1:length(f$x)){
  lines(c(f$x[i], f$x[i]), c(f$y[i] - 2 * sigma.n, f$y[i] + 2 * sigma.n), col = 'red', lwd = 2)
}