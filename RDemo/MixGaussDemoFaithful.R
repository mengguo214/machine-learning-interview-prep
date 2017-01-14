# Description: 
#   Illustration of the EM algorithm (1-D) for Gaussian Mixture Model (GMM)
#   using the Old Faithful dataset.
#
# Details:
#    K: the number of clusters is set to 2.
#
# License: GNU v3
#
# Date: October 19, 2016
#
# Authors:
#    Originally called EM Mixture MV.R.
#    Rcode was written by Michael Clark.
#    Rcode was modified by Jiali Lin.
#    Depts. of Statistics, Virginia Tech,
#    Hutcheson Hall, 403K, Blacksburg, VA 24061 
#
# Sources:
#    https://github.com/m-clark/Miscellaneous-R-Code/tree/master/ModelFitting/EM%20Examples
# 
# See Also:
#     http://www.geyserstudy.org/geyser.aspx?pGeyserNo=OLDFAITHFUL

rm(list = ls())

GaussmixEM <- function(params, X, clusters = 2, tol = .00001, maxits = 100, showits = T){     
  # EM for gaussian mixture
  #
  # Usage:
  #    GaussmixEM(params, X)
  # 
  # Input:
  #   params: initial parameters (means, covariances, cluster probability)
  #   X: data
  #   clusters: number of clusters desired
  #   tol: tolerance
  #   maxits: maximum iterations
  #   showits: whether to show iterations
  #  
  # Output:
  #   probs: probability
  #   mu: mean
  #   var: variance 
  #   resp: responsibilities
  #   cluster: cluster assignment 
  
  # Starting points
  N <- nrow(X)
  nams <- names(params)
  mu <- params$mu
  var <- params$var
  probs <- params$probs
  
  # Other initializations
  # Initialize cluster 'responsibilities', i.e. probability of cluster membership for each observation i
  ri <- matrix(0, ncol = clusters, nrow = N)  
  it <- 0
  converged <- FALSE
  
  # Show iterations
  if (showits)                                  
    cat(paste("Iterations of EM:", "\n"))
  
  while ((!converged) & (it < maxits)) { 
    probsOld <- probs
    muOld <- mu
    varOld <- var
    riOld <- ri
    
    # E step
    # ------------
    # Compute responsibilities
    for (k in 1:clusters){
      ri[, k] <- probs[k] * dnorm(X, mu[k], sd = sqrt(var[k]), log = F)
    }
    ri <- ri / rowSums(ri)
    
    # M step
    # ------------
    # rk is the weighted average cluster membership size
    rk <- colSums(ri)                           
    probs <- rk/N
    mu <- (t(X) %*% ri) / rk                      
    var <- (t(X^2) %*% ri) / rk - mu^2
    
    parmlistold <- rbind(probsOld, muOld, varOld)
    parmlistcurrent <- rbind(probs, mu, var)
    it <- it + 1
    
    # if showits true, & it =1 or divisible by 5 print message
    if (showits & it == 1 | it%%5 == 0){
      cat(it, "\n")
    }        
    converged <- max(abs(parmlistold - parmlistcurrent)) < tol
  }
  
  # create cluster membership 
  clust <- which(round(ri) == 1, arr.ind = T)    
  
  # order accoring to row rather than cluster
  clust <- clust[order(clust[, 1]), 2]           
  out <- list(probs = probs, mu = mu, var = var, resp = ri, cluster = clust)
} 

# faithful data set
data(faithful)

# starting parameters; requires mean, variance and class probabilitiy
params1 <- list(mu = c(2, 5), var = c(1, 1), probs = c(.5, .5)) 
params2 <- list(mu = c(50, 90), var = c(1, 15), probs = c(.5, .5))  

X1 <- matrix(faithful[, 1])
X2 <- matrix(faithful[, 2])

test1 <- GaussmixEM(params1, X = X1, tol = 1e-8)
test2 = GaussmixEM(params2, X = X2, tol = 1e-8)

library(ggplot2)
qplot(x = eruptions, y = waiting, data = faithful) 

# clustering by "eruptions"
ggplot(aes(x = eruptions, y = waiting), data = faithful) +
  geom_point(aes(color = factor(test1$cluster)))  +
  theme(legend.title = element_blank())

# clustering by "waiting"
ggplot(aes(x = eruptions, y = waiting), data = faithful) +
  geom_point(aes(color = factor(test2$cluster)))  +
  theme(legend.title = element_blank())

