# Description: 
#   This script illustrates the usage of the function 'GmmVbem.r'
#
# License: GNU v3
#
# Date: October 19, 2016
#
# Authors:
#    Originally called gmmVBEM.m.
#    Rcode was written by Emtiyaz, CS, UBC.
#    Rcode was modified by Jiali Lin.
#    Depts. of Statistics, Virginia Tech,
#    Hutcheson Hall, 403K, Blacksburg, VA 24061 
#
# References:
#   http://www.cs.ubc.ca/~murphyk/Software/VBEMGMM/index.html
#   http://stackoverflow.com/questions/28317203/r-superimposing-bivariate-normal-density-ellipses-on-scatter-plot

rm(list = ls())

if(!require("mixtools")) { install.packages("mixtools");  require("mixtools") }
library(matlab)

GmmVbem = function(x, mix, PriorPar, options){
  # Variational Bayes EM algorithm for Gaussian Mixture Model 
  
  # INPUT:
  #   D: the dimension;
  #   N: the number of Data points;
  #   x:   training data of size DxN ;
  #   mix: gmm model initialize with netlab's GmmVBEM function;
  #   PriorPar: structure containing priors;
  #   options: options for maxIter, threshold etc. etc.
  
  # initialize variables
  D = dim(x)[1];
  N = dim(x)[2];
  K = mix$ncentres;
  eps = 2.2204e-16
  likIncr = options$threshold + eps;
  L = rep(likIncr, options$maxIter);
  logLambdaTilde = matrix(0,1,K);
  E = matrix(0,nrow=N,ncol=K);
  trSW = matrix(0,1,K);
  xbarWxbar = matrix(0,1,K);
  mWm = matrix(0,1,K);
  trW0invW = matrix(0,1,K);
  
  # priors
  # prior over the mixing coefficients - CB 10.39
  alpha0 = PriorPar$alpha;
  # prior over gaussian means -  CB 10. 40
  m0 = PriorPar$mu;
  # prior over the gaussian variance - CB 10.40
  beta0 = PriorPar$beta;
  # wishart prior variables - CB 10.40
  W0 = PriorPar$W;
  W0inv = solve(W0);
  v0 = PriorPar$v;
  
  # Use 'responsibilities' from initialization to set sufficient statistics - 
  # CB 10.51-10.53.
  Nk = N*t(mix$priors);
  xbar = t(mix$centres);
  S = mix$covars;
  
  # Use above sufficient statistics for M step update equations - CB 10.58, 
  # 10.60-10.63. 
  alpha = alpha0 + Nk;
  beta = beta0 + Nk;
  v = v0 + Nk;
  m = ((beta0*m0)%*%matrix(1,1,K) + (matrix(1,D,1)%*%t(Nk))*xbar)/(matrix(1,D,1)%*%t(beta));
  W = array(0, c(D,D,K));
  
  for (k in 1:K) {
    mult1 = beta0*Nk[k]/(beta0 + Nk[k]);
    diff3 = xbar[,k] - m0;
    W[ , , k] = solve(W0inv + Nk[k]*S[ , ,k]  + mult1*diff3%*%t(diff3));
  }
  
  # Main loop of algorithm
  for (iter in 1:options$maxIter){
    # Calculate r - CB 10.64 - 10.66
    psiAlphaHat = psigamma(sum(alpha),0);
    logPiTilde = psigamma(alpha,0) - psiAlphaHat;
    const = D*log(2);
    
    for (k in 1:K){
      t1 = psigamma(0.5*repmat(v[k]+1,D,1) - 0.5*matrix(1:D), 0);
      logLambdaTilde[k] = sum(t1) + const + log(det(W[ , , k]));
      
      for (n in 1:N){
        # Calculate E
        diff = x[,n] - m[,k];
        E[n,k] = D/beta[k] + v[k]*t(diff)%*%W[ , , k]%*%diff;
      } # end of n
    } # end of k
    
    # Calculate rho - CB 10.45 - 10.46
    logRho = repmat(matrix(t(logPiTilde) + 0.5*logLambdaTilde),c(N,1)) - 0.5*E;
    logSumRho = apply(logRho, 1, function(x) log(sum(exp(x)))); #?
    logr = logRho - repmat(logSumRho, 1,K);
    r = exp(logr);
    
    # compute N(k) - CB 10.51
    Nk = exp(apply(logr, 2, function(x) log(sum(exp(x))))) # not matrix yet 15*1
    # add a non-zero term for the components with zero responsibilities
    Nk = Nk + 1e-10; 
    # compute xbar(k), S(k) - CB 10.52 - 10.53
    for (k in 1:K){
      xbar[ ,k] = rowSums(repmat(matrix(t(r[,k])),c(D,1))*x)/Nk[k];
      diff1 = x - repmat(xbar[,k],1,N);
      diff2 = repmat(matrix(t(r[,k])),c(D,1))*diff1;
      S[, , k] = (diff2%*%t(diff1))/Nk[k];
    }
    
    # compute Lower bound (refer to Bishop for these terms) - CB 10.71 - 10.77
    # C(alpha0)
    logCalpha0 = lgamma(K*alpha0) - K*lgamma(alpha0);
    # B(lambda0)
    logB0 = (v0/2)*log(det(W0inv)) - (v0*D/2)*log(2) -
      (D*(D-1)/4)*log(pi) - sum(lgamma(0.5*(v0+1 -1:D)));
    # log(C(alpha))
    logCalpha = lgamma(sum(alpha)) - sum(lgamma(alpha));
    # Various other parameters for different terms
    H =0;
    for (k in 1:K){
      # sum(H(q(Lamba(k))))
      logBk = -(v[k]/2)*log(det(W[ , , k])) - (v[k]*D/2)*log(2)
      - (D*(D-1)/4)*log(pi) - sum(lgamma(0.5*(v[k] + 1 - 1:D)));
      H = H -logBk - 0.5*(v[k] -D-1)*logLambdaTilde[k] + 0.5*v[k]*D;
      # for Lt1 - third term
      trSW[k] = sum(diag(v[k]*S[ , , k]%*%W[ , , k]));
      diff = xbar[,k] - m[,k] ;
      xbarWxbar[k] = t(diff)%*%W[ , , k]%*%diff;
      # for Lt4 - Fourth term
      diff = m[,k] - m0;
      mWm[k] = t(diff)%*%W[ , , k]%*%diff;
      trW0invW[k] = sum(diag(W0inv%*%W[ , , k]));
    }
    
    #Bishop's Lower Bound
    Lt1 = 0.5*sum(Nk*(t(logLambdaTilde) - D/beta - t(trSW) - v*t(xbarWxbar) - D*log(2*pi)));
    Lt2 = sum(Nk*logPiTilde)
    Lt3 = logCalpha0 + (alpha0 -1)*sum(logPiTilde);
    Lt41 = 0.5*sum(D*log(beta0/(2*pi)) + t(logLambdaTilde) - D*beta0/beta - beta0*v*t(mWm));
    Lt42 = K*logB0 + 0.5*(v0-D-1)*rowSums(logLambdaTilde) - 0.5*sum(v*t(trW0invW));
    Lt4 = Lt41+Lt42;
    Lt5 = sum(sum(r*logr));
    Lt6 = sum((alpha - 1)*logPiTilde) + logCalpha;
    Lt7 = 0.5*sum(t(logLambdaTilde) + D*log(beta/(2*pi))) - 0.5*D*K - H;
    
    #Bishop's Lower Bound 
    L[iter] = Lt1 + Lt2 + Lt3 + Lt4 - Lt5 - Lt6 - Lt7; # where is L?
    
    # warning  if lower bound decreses
    if (iter>2 && L[iter]<L[iter-1]) {
      cat("Lower bound decreased by =", L[iter]-L[iter-1], "\n" )
    }
    
    # Begin M step
    # compute new parameters - CB 10.58 - 10.63
    alpha = alpha0 + Nk;
    beta = beta0 + Nk;
    v = v0 + Nk;
    m = (repmat(beta0*m0,1,K) + repmat(matrix(t(Nk)),c(D,1))*xbar)/repmat(matrix(t(beta)),c(D,1));
    for (k in 1:K) {
      mult1 = beta0*Nk[k]/(beta0 + Nk[k]);
      diff3 = xbar[,k] - m0;
      W[, , k]  = solve(W0inv + Nk[k]*S[, , k]  + mult1*diff3%*%t(diff3));
    }
    
    #PLOT 
    if (options$displayIter){
      cat("iter =", iter, "\n" )
    }
    
    if (options$displayFig){
      dev.off()
      plot(x[1,], x[2,], xlab = "eruptions", ylab = "waiting")
      for (i in 1:K)  ellipse(m[,i],solve(W[,,i])/(v[i]-D-1),col = i)
      Sys.sleep(0.1)
    }
    
    if (iter>1){
      likIncr = abs((L[iter]-L[iter-1])/L[iter-1]);
    }
    
    if (likIncr < options$threshold) break
  } #end iter
  
  return(list(L = L, center = m))
}

# Example of GMM VBEM clustering by "faithful" data
# ------------

# read data
X = faithful
plot(X)

# standardize the data
X = t(X)
X = X - repmat(apply(X,1,mean),1,ncol(X));
X = X/repmat(apply(X,1,var),1,ncol(X));
dim = dim(X)[1]
N = dim(X)[2]

# initialize with EM algorithm 
ncentres = 15;
mix = list();
mix$ncentres = ncentres;
mix$priors =  matrix(1,1,mix$ncentres)/ mix$ncentres;
mix$centres = matrix(rnorm(mix$ncentres*dim,mean=0,sd=1), 
                     mix$ncentres, dim);
mix$covars = repmat(diag(dim), c(1, 1, mix$ncentres));

# intialize the priors
PriorPar = list();
PriorPar$alpha = .001;
PriorPar$mu = matrix(0,dim,1);
PriorPar$beta = 1;
PriorPar$W = 200*diag(dim);
PriorPar$v = 20;

# set the options for VBEM
options = list();
options$maxIter = 100;
options$threshold = 1e-5;
options$displayFig = TRUE;
options$displayIter = TRUE;

# Call the function
res = GmmVbem(X, mix, PriorPar, options)

# plot lower bound
plot(res$L, lty = 2, type ="l", 
     ylab = "lower bound on log marginal likelihood",
     main = "variational Bayes objective for GMM on old faithful data")
