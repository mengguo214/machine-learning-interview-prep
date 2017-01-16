# Description: 
#   Variational bayes for univariate Gaussian
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

# Make data with mean 0 and unit variance
set.seed(12);
N = 10; D = 1;
data = runif(N);
data = (data - mean(data)) / sd(data);
m = mean(data);
s = sum(data);
data = matrix(data);
sSq = sum(data*data);
xbar = m;
sigma2Hat = sSq/N - xbar^2;

# The true posterior model using normal gamma prior
# hyper parameters
a0 = 0; b0 = 0; mu0 = 0; kappa0 = 0;
truePost = list()
truePost$mu = mean(data);
truePost$kappa = N;
truePost$alpha = N/2;
truePost$beta = 1/2*sum((data - m)^2);

# Initialize VB to fit univariate gaussian
# initial guess
aN = 2.5; bN = 1;
muN = 0.5; kappaN = 5;

# plot target distribution
mu <- seq(-1.0 , 1.0 , by = 0.1)
lambda <- seq(0, 2, by = 0.1)

TrueNormalGammaPdf <- function(mu, lambda){
  muprior = truePost$mu; kappa = truePost$kappa;  
  alpha = truePost$alpha; beta = truePost$beta;
  C = (beta^alpha * sqrt(kappa)) / (gamma(alpha) * sqrt(2 * pi))
  p = C * (lambda^(alpha-1/2)) * (exp(-beta * lambda)) *  
    (exp(-kappa/2* (lambda * (mu - muprior)^2)))
  return(p)
}

# Posterior
vbPost <- function(mu, lambda){
  C = (bN^aN * sqrt(kappaN)) / (gamma(aN) * sqrt(2 * pi))
  p = C * (lambda^(aN-1/2)) * (exp(-bN * lambda)) *  
    (exp(-kappaN/2* (lambda * (mu - muN)^2)))
  return(p)
}

dens <- outer(mu, lambda, TrueNormalGammaPdf)

par(mfrow = c(2, 2))
contour(mu, lambda, dens, 
        xlab = 'mu', ylab = 'lambda', col = 'green', drawlabels = F)

iter = 1;
Lq = 0; 
maxIter = 100;
converged = FALSE;
tol = 1e-5;
Lbound = rep(NA, maxIter); Lbound[1] = Lq;

while ((!converged) & (iter < maxIter)) {
  LqOld = Lq;
  
  # update q(mu)
  elambda = aN / bN;
  muN = (kappa0 * mu0 + N * m) / (kappa0 + N);
  kappaN = (kappa0 + N) * elambda;
  
  if (iter == 1){
    dens <- outer(mu, lambda, TrueNormalGammaPdf)
    contour(mu, lambda, dens, 
            xlab = 'mu', ylab = 'lambda', col = 'green', drawlabels = F)
    dens <- outer(mu, lambda, vbPost)
    contour(mu, lambda, dens, 
            xlab = 'mu', ylab = 'lambda', col = 'red', add = T, drawlabels = F)
  }
  
  # update q(lambda)
  emu = muN;
  emuSquare = 1/kappaN + muN^2;
  aN = a0 + (N+1)/2;
  bN = b0 + 1/2 * ((sSq + kappa0 * mu0^2) - 2 * emu *(s + kappa0 * mu0) + emuSquare*(kappa0+N));
  if (iter == 1){
    dens <- outer(mu, lambda, TrueNormalGammaPdf)
    contour(mu, lambda, dens, 
            xlab = 'mu', ylab = 'lambda', col = 'green', drawlabels = F)
    dens <- outer(mu, lambda, vbPost)
    contour(mu, lambda, dens, 
            xlab = 'mu', ylab = 'lambda', col = 'red', add = T, drawlabels = F)
  }
  
  Lq = 1/2 * log(1/kappaN) + log(gamma(aN)) - aN * log(bN);
  Lbound[iter] = Lq;
  if (iter > 1){
    if (Lq - LqOld < -tol){
      cat('Lq did not increase, iter =', iter, "\n");
    } else if (abs(Lq - LqOld) < tol){
      converged = TRUE;
    }
  }
  
  iter = iter + 1;
  
}

cat('Total # of iterations:', iter, "\n");
dens <- outer(mu, lambda, TrueNormalGammaPdf)
contour(mu, lambda, dens, 
        xlab = 'mu', ylab = 'lambda', col = 'green', drawlabels = F)
dens <- outer(mu, lambda, vbPost)
contour(mu, lambda, dens, 
        xlab = 'mu', ylab = 'lambda', col = 'red', add = T, drawlabels = F)