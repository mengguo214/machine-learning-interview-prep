# Description: 
#   This script shows simulation of multiple linear regression 
#   with variable Selection, using Gibbs sampler. 
#
# Details:
#    Need to know the full conditional distribution of parameters.
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
#    http://www.apps.stat.vt.edu/zhu/teaching/2016/6474/6474_2016.htm

rm(list = ls())
library(mvtnorm)
library(pscl)
library(latex2exp)
library(ggplot2)

BvsGibbs <-  function(covariate, Y, nMc, v1, v0, w, gamma0, sigma0, a, b){
  
  # Bayesian Variables Selection:
  # Continuous Gaussian Responses Gibbs Sampling.
  #
  # Input: 
  #   covariate: n by p covariate matrix (does not include the column of ones). 
  #   Y: continuous response, n by 1. 
  #   nMc: number of MCMC iterations.
  #   v1: p by 1 vector, the large s.d. of beta when gamma_j=1.
  #   v0: p by 1 vector, the small s.d. of beta when gamma_j=0.
  #   w: the probablity that gamma_j=1 for each variable (=expected num. of variables selected/p).
  #   gamma0: initial value of gamma, contains zeros and ones.
  #   beta0: initial value of beta.
  #   sigma0: initial value of sigma^2. 
  #   a, b: prior parameters for sigma^2 (Inverse Gamma parameters).
  #
  # Output:
  #   sample.beta: All posterior samples of beta.
  #   sample.gamma: all posterior samples of gamma.
  #   sample.sigma2: All posterior samples of sigma2.
  #   logodds: The log (conditional) posterior odds used when updating gamma.
  
  n <-  dim(covariate)[1]        
  p <-  dim(covariate)[2]
  X <-  cbind(1,covariate) 
  
  library(mvtnorm)	
  save.beta <-  matrix(NA, nMc, p + 1) 
  save.gamma <-  matrix(NA, nMc, p)    
  save.odds <-  matrix(NA, nMc, p)   
  save.sigma2 <-  rep(NA, nMc)
  
  gamma <-  gamma0
  sigma2 <-  sigma0
  
  for (i in 1:nMc){
    
    # Update beta.
    invSIG <-  diag(c(1 / 10^6, gamma * (1 / v1^2) + (1 - gamma) * (1 / v0^2))) #assume R=I. 
    K <-  t(X) %*% X / sigma2 + invSIG
    invK <-  solve(K)
    M <-  (t(X) %*% Y) / sigma2
    beta <-  t(rmvnorm(n = 1, mean = invK%*%M, sigma = invK)) 	
    save.beta[i, ] <-  beta
    
    # Update gamma from (conditional) posterior odds.
    gamma <-  rep(NA, p) 
    log.odds <-  -log(v1) - beta[2:(p + 1)]^2 / (2 * v1^2) + log(v0) + 
      beta[2:(p + 1)]^2 / (2 * v0^2) + log(w / (1 - w))
    
    save.odds[i, ] <-  log.odds
    for (j in 1:p){
      if (log.odds[j] > 0)  gamma[j] <-  1
      else {        
        gamma[j] <-  rbinom(1, size = 1, prob = exp(log.odds[j])/(exp(log.odds[j]) + 1))        
      }
    }
    save.gamma[i,] <-  gamma
    
    # Update sigma2. 
    a.tilde <-  a + n/2
    b.tilde <-  sum((Y - X %*% beta)^2) / 2 + b
    sigma2 <-  rgamma(1, a.tilde, b.tilde)
    save.sigma2[i] <-  sigma2
  }
  return(list(sample.beta = save.beta, sample.gamma = save.gamma, 
              sample.sigma2 = save.sigma2, logodds = save.odds))
}

# Generate data.
set.seed(50)
sigma2 <-  4
p <-  5 
n <-  100 
true.gamma <-  c(0, 1, 0, 0, 1)  
true.beta <-  c(10, 0, 8, 0, 0, 6)   
covariate <-  matrix(rnorm(n*p, 0, 1), ncol = p, byrow = F) 
Y <-  cbind(1, covariate) %*% true.beta + rnorm(n, 0, sqrt(sigma2)) 

# Rough estimate without initial values
X0 <-  cbind(1, covariate)
beta.hat <-  solve(t(X0) %*% X0) %*% t(X0) %*% Y
sig2.hat <-  sum((Y - X0 %*% beta.hat)^2) / n

# Set inital values and parameters for MCMC
nMc <-  10000
a <-  (sig2.hat)^2/10 + 2 # determine (a,b) s.t. E(sig2)=sig2.hat, var(Sig2)=10
b <-  sig2.hat * (a - 1)
v0 <-  rep(10e-3, p)
v1 <-  rep(sqrt(max(beta.hat) * 100), p)
w <-  1 - sum(abs(beta.hat) < 4) / p
gamma0 <-  rbinom(p, size = 1, prob = w)
sigma0 <-  sig2.hat

# Tic
ptm <-  proc.time()

# Run MCMC
res <-  BvsGibbs(covariate, Y, nMc, v1, v0, w, gamma0, sigma0, a, b)

# Tok
Time.used <-  proc.time() - ptm
Time.used[3]/60

# Get posterior mean
burnin <-  3000
beta.sample <-  res$sample.beta[(burnin + 1):nMc, ]
beta.mean <-  colMeans(beta.sample)
par(mfrow <-  c(1, 1))

# Posterior inference of beta.
qplot(x = 1:6, y = true.beta, colour = "Truth", geom = "point") + 
  geom_point(aes(y = beta.mean, colour = "Estimation")) + 
  ylab(TeX("$\\gamma_j$ vs $\\hat{\\gamma_j}$")) + 
  theme(legend.title = element_blank())

# Post prob. of gamma_j=1
gamma.sample <-  res$sample.gamma[(burnin + 1):nMc, ]
post.gamma <-  colMeans(gamma.sample)
qplot(x = 1:5, y = true.gamma, colour = "Truth", geom = "point") + 
  geom_point(aes(y = post.gamma, colour = "Estimation")) + 
  coord_cartesian(ylim = c(0, 5)) + 
  coord_cartesian(xlim = c(1, 5)) +
  ylab(TeX("$\\gamma_j$ vs $\\hat{\\gamma_j}$")) + 
  theme(legend.title = element_blank())

# Post. mean of sigma2
sigma2.sample <-  res$sample.sigma2[(burnin + 1):nMc]
sig2.mean <-  mean(sigma2.sample)  
qplot(sigma2.sample, geom = "histogram") + 
  geom_vline(xintercept = sig2.mean, color = "blue") +
  xlab(TeX("$\\hat{\\sigma^2}$ sample"))  

