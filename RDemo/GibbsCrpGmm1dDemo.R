# Description: 
#   This script shows univariate Gaussian DP mixture model
#
# Details:
#   Likelihood: x_i|mu_i,T_i ~ N(mu, 1/T)
#   Random effect: (mu, T) ~ G
#   Prior: G ~ DP(alpha, G0), where G0 ~ NormalGamma(mu0, kappa0, tau0, beta0)
#   
#   Posterior: mu, T|X ~ NormalGamma(muN, kappaN, tauN, betaN), where
#     xm: mean(x)
#     muN: (kappa0 * mu0 + N * xm) / kappaN
#     kappaN: kappa0 + N
#     tauN: tau0 + N / 2
#     betaN: beta0 + 0.5 * sum((x - xm)^2) + (kappa0 * N *(xm - mu0)^2) / (2 * kappaN)
#
#   Marginal Distributions:
#   By construction, the marginal distribution over T is a gamma distribution.
#   The conditional distribution over X given T is a Gaussian distribution. 
#   The marginal distribution over X is a three-parameter non-standardized Student's t-distribution

# License: GNU v3
#
# Date: October 19, 2016
#
# Authors:
#    Originally called dpmm.R.
#    Rcode was written by John Myles White.
#    Rcode was modified by Jiali Lin.
#    Depts. of Statistics, Virginia Tech,
#    Hutcheson Hall, 403K, Blacksburg, VA 24061 
#
# Source: 
#   https://github.com/johnmyleswhite/bayesian_nonparametrics/blob/master/code/dpmm/dpmm.R
# 
# Reference:
#   https://en.wikipedia.org/wiki/Normal-gamma_distribution
#   https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf

rm(list = ls())

library(ggplot2)

t_logpdf <- function(x, mu, v, df)
{
  # Generalized t-Distribution log PDF
  g <- lgamma(0.5 * df + 0.5) - lgamma(0.5 * df) - log(sqrt(df * pi * v))
  return(g - 0.5 * (df + 1) * log(1 + (1 / df) * ((x - mu)^2) / v))
}

normalize_logprob <- function(log.prob)
{
  # Normalize a log probability vector without numerical underflow/overflow
  
  # Compute the log of the normalizing constant.
  g <- log(sum(exp(log.prob - max(log.prob)))) + max(log.prob)
  
  # Find probabilities by exponentiating normalized log probabilities.
  return(exp(log.prob - g))
}

dpmm_loglik <- function(xn, x, tau0, beta0, mu0, kappa0)
{
  # Calculate the log likelihood for the Gaussian DP mixture model
  # Mean and variance parameters marginalized under normal-gamma prior
  # This corresponds to a generalized Student's t-distribution
  # NB: some constant terms are ignored
  
  N <- length(x)
  kappaN <- kappa0 + N
  tauN <- tau0 + N / 2
  
  if (N > 0)
  {
    # If there previous observations, use the current posterior.
    xm <- mean(x)
    betaN <- beta0 + 0.5 * sum((x - xm)^2) + (kappa0 * N *(xm - mu0)^2) / (2 * kappaN)
    muN <- (kappa0 * mu0 + N * xm) / kappaN
  }
  else
  {
    # If there are no previous observations, revert to the prior.
    betaN <- beta0
    muN <- mu0
  }
  
  return(t_logpdf(xn, muN, betaN * (kappaN + 1) / (tauN * kappaN), 2 * kappaN))
}


dpmm_gibbs <- function(x, alpha, tau0, beta0, mu0, kappa0, nIter)
{
  # Gibbs sampling for the Gaussian Dirichlet process mixture model
  #
  # Inputs:
  #  x - data (vector of length N)
  #  alpha - DP concentration parameter
  #  tau0 - normal-gamma prior shape
  #  beta0 - normal-gamma prior rate (inverse scale)
  #  kappa0 - normal-gamma prior precision scaling parameter
  #  nIter - number of Gibbs iterations
  #
  # Outputs:
  #  C, a N x nIter matrix of cluster assignments
  
  N <- length(x)
  p <- rep(1, 1, N)
  C <- matrix(data = NA, nrow = N, ncol = nIter); c <-rep(1, 1, N)
  m <- rep(0, 1, N); m[1] <- N
  logpost <- rep(1, 1, nIter); logprior <- rep(1, 1, N); loglik <- rep(1, 1, N)
  ix <- 1:N
  
  for (i in 1:nIter)
  {
    print(paste("Iteration", i))
    
    for (n in 1:N)
    {
      #all customers except n
      cn <- c[1:N!=n]
      
      #count cluster assignments
      for (k in 1:N)
      {
        m[k] <- sum(cn == k)
      }
      
      if (all(m > 0))
      {
        #active dishes
        K.active <- ix[m > 0]
      }
      else
      {
        #active dishes + 1 new dish
        K.active <- c(ix[m > 0], min(ix[m == 0]))
      }
      for (k in K.active)
      {
        if (m[k] > 0)
        {
          #prior for old dish
          logprior[k] <- log(m[k])
        }
        else
        {
          #prior for new dish
          logprior[k] <- log(alpha)
        }
        #calculate log likelihood
        loglik[k] <- dpmm_loglik(x[n], x[c==k], tau0, beta0, mu0, kappa0)
      }
      
      #posterior
      post <- normalize_logprob(logprior[K.active] + loglik[K.active])
      
      #update cluster assignment
      c[n] <- sample(K.active, 1, rep = TRUE, prob = post)
    }
    
    C[,i] = c
  }
  C
}

set.seed(1)
x <- c(rnorm(100, 100, 8), rnorm(50, 200, 25), rnorm(50, 25, 1))
labels <- c(rep("A", 100), rep("B", 50), rep("C", 50))

df <- data.frame(X = x, Label = labels)

ggplot(df, aes(x = X)) +
  geom_histogram(binwidth = 3)

ggplot(df, aes(x = X, fill = Label)) +
  geom_histogram(binwidth = 3) +
  labs(legend.position = "none") +
  labs(title = "Ground Truth Mixture Model")
#ggsave("dpmm_ground_truth.pdf", height = 7, width = 7)

nIter <- 100
results <- dpmm_gibbs(x, 0.01, 0.1, 0.1, 0, 0.1, nIter)
results[, nIter]

ggplot(df, aes(x = X, fill = factor(results[, nIter]))) +
  geom_histogram(binwidth = 3) +
  labs(legend.position = "none") +
  labs(title = "dp-MM with alpha = 0.01")
#ggsave("dpmm_0.01.pdf", height = 7, width = 7)

nIter <- 100
results <- dpmm_gibbs(x, 0.5, 0.1, 0.1, 0, 0.1, nIter)
results[, nIter]

ggplot(df, aes(x = X, fill = factor(results[, nIter]))) +
  geom_histogram(binwidth = 3) +
  labs(legend.position = "none") +
  labs(title = "dp-MM with alpha = 0.5")
#ggsave("dpmm_0.5.pdf", height = 7, width = 7)

nIter <- 100
results <- dpmm_gibbs(x, 2.5, 0.1, 0.1, 0, 0.1, nIter)
results[, nIter]

ggplot(df, aes(x = X, fill = factor(results[, nIter]))) +
  geom_histogram(binwidth = 3) +
  labs(legend.position = "none") +
  labs(title = "dp-MM with alpha = 2.5")
#ggsave("dpmm_2.5.pdf", height = 7, width = 7)

nIter <- 100
results <- dpmm_gibbs(x, 12.5, 0.1, 0.1, 0, 0.1, nIter)
results[, nIter]

ggplot(df, aes(x = X, fill = factor(results[, nIter]))) +
  geom_histogram(binwidth = 3) +
  labs(legend.position = "none") +
  labs(title = "dp-MM with alpha = 12.5")
#ggsave("dpmm_12.5.pdf", height = 7, width = 7)

nIter <- 100
results <- dpmm_gibbs(x, 100.0, 0.1, 0.1, 0, 0.1, nIter)
results[, nIter]

ggplot(df, aes(x = X, fill = factor(results[, nIter]))) +
  geom_histogram(binwidth = 3) +
  labs(legend.position = "none") +
  labs(title = "dp-MM with alpha = 100.0")
#ggsave("dpmm_100.0.pdf", height = 7, width = 7)

