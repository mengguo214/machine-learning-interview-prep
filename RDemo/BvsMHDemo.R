# Description: 
#   This script shows simulation of multiple linear regression 
#   with variable selection, using MH sampler.
#
# Details:
#    For MH sampling scheme, we use switch and swap proposal (Brown, Vannucci, Fearn (1998), JRSS B).
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

BvsMH <- function(covariate,Y,nMc,v1,v0, w, gamma0, a, b, phi){

  # Bayesian Variables Selection:
  # Continuous Gaussian Responses MH Sampling. 
  #
  # Input: 
  #   covariate: n by p covariate matrix (does not include the column of ones). 
  #   Y: continuous response, n by 1. 
  #   nMc: number of MCMC iterations
  #   v1: p by 1 vector, the large s.d. of beta when gamma_j=1
  #   v0: p by 1 vector, the small s.d. of beta when gamma_j=0
  #   w: the probablity that gamma_j=1 for each variable (=expected num. of variables selected/p)
  #   a, b: prior parameters for sigma^2. (Inverse Gamma parameters)
  #   gamma0:initial value of tao, contains zeros and ones
  #
  # Output:
  #   sample.gamma: all posterior samples of gamma

  
  n <- dim(covariate)[1]        
  p <- dim(covariate)[2]
  X <- cbind(1,covariate) 
  
  save.gamma <- matrix(NA, nMc, p)   
  save.ratio <- rep(NA, nMc)    
  gamma <- gamma0
  phi <- 0.5
  accept.count <- 0
  
  for (i in 1:nMc){
    
    #Update gamma using stochastic search.
    
    # propose a new gamma
    gamma.n <- gamma
    if (sum(gamma) == p)	 gamma.n[sample(1:p, 1)] <- 0 
    if (sum(gamma) == 0) 	 gamma.n[sample(1:p, 1)] <- 1
    if ((sum(gamma) < p) & (sum(gamma) > 0)){
      if (runif(1) > phi) {
        chid <- sample(1:p, 1)
        gamma.n[chid] <- 1 - gamma[chid] # switch with prob. phi. 
        
      } else{ 
        l0 <- which(gamma == 0)
        l1 <- which(gamma == 1)
        gamma.n[l0[sample(1:length(l0), 1)]] <- 1
        gamma.n[l1[sample(1:length(l1), 1)]] <- 0 
      }
    } 
    
    Sig.n <- diag(c(10^6, gamma.n * v1^2 + (1 - gamma.n) * v0^2)) # assume that R=I.
    invSig.n <- diag(c(1/10^6, gamma.n * (1 / v1^2) + (1 - gamma.n) * (1 / v0^2)))
    A.new <- b + (sum(Y^2) - t(Y) %*% X %*% solve(t(X) %*% X + invSig.n) %*% t(X) %*% Y) / 2
    
    Sig.old <- diag(c(10^6, gamma * v1^2 + (1 - gamma) * v0^2))
    invSig.old <- diag(c(1 / 10^6, gamma * (1 / v1^2) + (1 - gamma)*(1 / v0^2)))
    A.old <- b + (sum(Y^2) - t(Y) %*% X %*% solve(t(X) %*% X + invSig.old) %*% t(X) %*% Y) / 2
    
    log.R <- -log(det(Sig.n %*% (t(X) %*% X) + diag(p + 1))) / 2 - (n / 2 + a) * log(A.new) + 
      log(det(Sig.old %*% t(X) %*% X + diag(p + 1))) / 2 + (n / 2 + a) * log(A.old) + 
      (sum(gamma.n) - sum(gamma)) * log(w / (1 - w))    
    
    save.ratio[i] <- log.R
    if (log(runif(1)) < log.R) {
      gamma <- gamma.n
      accept.count <- accept.count + 1} 
    
    save.gamma[i,] <- gamma
    
  }
  
  return(list(sample.gamma = save.gamma, 
              post.ratio = save.ratio, 
              accept.count = accept.count))
}

# Generate data.
set.seed(50)
sigma2 <- 4
p <- 5 
n <- 100 
true.gamma <- c(0, 1, 0, 0, 1)  
true.beta <- c(10, 0, 8, 0, 0,  6)   
covariate <- matrix(rnorm(n*p, 0, 1), ncol = p, byrow = F) 
Y <- cbind(1, covariate) %*% true.beta + rnorm(n, 0, sqrt(sigma2)) 

# Rough estimate without initial values
X0 <- cbind(1, covariate)
beta.hat <- solve(t(X0) %*% X0) %*% t(X0) %*% Y
sig2.hat <- sum((Y - X0 %*% beta.hat)^2) / n

# Set inital values and parameters for MCMC
nMc <- 10000
a <- (sig2.hat)^2/10 + 2 # determine (a,b) so that E(sig2)=sig2.hat, var(Sig2)=10
b <- sig2.hat * (a - 1)
v0 <- rep(10e-6, p)
v1 <- rep(sqrt(max(beta.hat) * 100), p)/sqrt(sig2.hat)
w <- 1 - sum(abs(beta.hat) < 4) / p
gamma0 <- rep(0, p)

# Tic 
ptm <- proc.time()

# Run MCMC
res <- BvsMH(covariate, Y, nMc, v1, v0, w, gamma0, a, b, phi)

# Tok
Time.used <- proc.time() - ptm
Time.used[3]/60

# Posterior inference of gamma.
burnin <- 3000
accep.rate <- res$accept.count/nMc
gamma.sample <- res$sample.gamma[(burnin + 1):nMc, ]
post.gamma <- colMeans(gamma.sample)

qplot(x = 1:5, y = true.gamma, colour = "Truth", geom = "point") + 
  geom_point(aes(y = post.gamma, colour = "Estimation")) + 
  coord_cartesian(ylim = c(0, 5)) + 
  coord_cartesian(xlim = c(1, 5)) + 
  xlab("x") + 
  ylab("posterior gamma") + 
  ggtitle(TeX("marginal posterior prob of $\\gamma_j=1$"))

