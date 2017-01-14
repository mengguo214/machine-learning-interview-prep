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
# Source:
#   http://www.apps.stat.vt.edu/zhu/teaching/2016/6474/6474_2016.htm

# References:
#    https://pmtk3.googlecode.com.

rm(list = ls())

# generate data from true distribution
set.seed(12234)
mu.true <- 0
lambda.true <- 0.8
N <- 100
x <- rnorm(N , mu.true, 1/lambda.true)

# initial value
# lambda ~gamma(a0,b0)
# mu ~ N(mu0,1/(kappa0*lambda))
mu0 <- 0.5
kappa0 <- 2
a0 <- 8
b0 <- 4

# specify the tolerance e
e <- 0.01

EA <- vector()
kappa <- vector()
b <- vector()

# an and mun are always the same
an <- a0 + (N + 1)/2
mun <- (kappa0 * mu0 + sum(x)) / (kappa0 + N)

# Main loop of algorithm
i <- 1
# step1
kappa[i] <- (kappa0 + N) * (a0 / b0)
#q(mu) ~ N(mu[1],1/kappa[1])

# step 2
# EA[i] = N*(mun^2 +1/kappa[i])+sum(x^2)-2*sum(x)*mun+kappa0/kappa[i]
EA[i] <- N * (mun^2 + 1 / kappa[i]) + sum(x^2) - 2 * sum(x) * mun + kappa0 ** (1 / kappa[i] + (mun - mu0)^2)
b[i] <- 0.5 * EA[i] + b0
# q(lambda)~ gamma(a[1], b[1])

# i=2
repeat{
  # step1 
  i <- i + 1
  kappa[i] <- (kappa0 + N) * (an / b[i-1])  
  #q(mu) ~ N(mu[1],1/kappa[1])
  
  # step 2
  # EA[i] = N*(mun^2 +1/kappa[i])+sum(x^2)-2*sum(x)*mun+kappa0/kappa[i]
  EA[i] <- N * (mun^2 + 1 / kappa[i]) + sum(x^2) - 2* sum(x) * mun + kappa0 * (1 / kappa[i] + (mun - mu0)^2)
  b[i] <- 0.5 * EA[i] + b0
  # q(lambda)~ gamma(a[1], b[1])
  
  if(abs(kappa[i] - kappa[i-1]) <e & abs(b[i] - b[i-1]) < e) break()
  kappa.optimal <- kappa[i]
  b.optimal <- b[i]
}


# Contour plot
# true posterior density p(mu,lambda/x)
f.true <- Vectorize(function(mu , lambda){
  f <- (sqrt(lambda/2 / pi))^N *exp(-(lambda / 2) * sum((x - mu)^2)) *
    sqrt(kappa0 * lambda / 2 /pi) * exp(-(kappa0 * lambda / 2)*(mu - mu0)^2 - b0 * lambda) * lambda^(a0 - 1)
  #logf=N/2*log(lambda)-(lambda/2)*sum((x-mu)^2)+0.5*log(kappa0*lambda)-0.5*kappa0*lambda*(mu-mu0)^2+(a0-1)*log(lambda)-b0*lambda
  
  #return(logf)
  return(f)
})
mu <- seq(-1.0 , 1.0 , by = 0.1)
lambda <- seq(0.0, 2,by = 0.1)
dens <- outer(mu, lambda, f.true)
#logdens=outer(mu,lambda,f.true)
#dens= exp(logdens - max(logdens))
par(mfrow = c(2, 2))
#  Contours of the true posterior distribution f(mu,lambda/X) --green

contour1<-c(1e-84, 1e-80, 1e-79, 5e-78, 1e-76)


#contours2 <- c(1.4e-100,1e-80,1e-50)
#par(mfcol=c(1,1))
contour(mu, lambda, dens, levels = contour1, xlab = 'mu', ylab = 'lambda', col = 'green', drawlabels = F)

#(a) Contours of the initial factorized approximation q(mu)*q(lambda) --blue
fa  <- Vectorize(function(mu, lambda){
  #   f.a = sqrt(kappa0*lambda/2/pi)*exp(-kappa0*lambda*0.5*(mu-mu0)^2)*lambda^(a0-1)*exp(-b0*lambda)  
  f.a <- (kappa0 * 0.5 * (mu - mu0)^2 + b0)^(-a0 - 0.5) * lambda^(a0 - 1) * exp(-b0 * lambda) 
  
  #  return(logfa)
  return(f.a)
})

#logdens.a=outer(mu,lambda,fa)
#dens.a= exp(logdens.a - max(logdens.a))
dens.a <- outer(mu, lambda, fa)
contour2 < -c(8e-7, 1e-7, 2e-7, 3e-7, 3.4e-7)

contour(mu, lambda, dens.a, levels = contour2, xlab = 'mu', ylab = 'lambda', col = 'blue', add = T, drawlabels = F)
legend("topleft",'(a)', bty = 'n')


#(b) After re-estimating the factor q(mu)--blue
fb  <- Vectorize(function(mu, lambda){
  f.b <- sqrt(kappa[1] * 0.5 / pi) * exp(-kappa[1] / 2 * (mu - mun)^2) * lambda^(a0-1) * exp(-b0 * lambda) 
  #  logfb=-kappa[1]*0.5*(mu-mun)^2+(a0-1)*log(lambda)-b0*lambda
  
  #  return(logfb)
  return(f.b)
})

#logdens.b=outer(mu,lambda,fb)
#dens.b= exp(logdens.b - max(logdens.b))
dens.b <- outer(mu, lambda, fb)
contour3 <- c(1e-4, 0.01, 0.1, 0.2)
contour(mu, lambda, dens, levels = contour1, xlab = 'mu', ylab = 'lambda', col = 'green', drawlabels = F)

contour(mu, lambda,dens.b, levels = contour3, add = T, drawlabels = F,col = 'blue')
legend("topleft",'(b)',bty = 'n')

#(c) After re-estimating the factor q(lambda)--blue
fc <- Vectorize(function(mu, lambda){
  f.c <- sqrt(kappa[1] / 2 / pi) * exp(-kappa[1] / 2 * (mu - mun)^2) * lambda^(an - 1) * exp(-b[1] * lambda)
  # logfc=-kappa[1]*0.5*(mu-mun)^2+(an-1)*log(lambda)-b[1]*lambda
  
  #  return(logfc)
  return(f.c)
})

#logdens.c=outer(mu,lambda,fc)
#dens.c= exp(logdens.c - max(logdens.c))
dens.c <- outer(mu, lambda, fc)
contour4 <- c (1e-45, 5e-43, 1e-40, 1e-38, 1e-37, 1e-35)
contour(mu, lambda, dens, levels = contour1, xlab = 'mu', ylab = 'lambda', col = 'green', drawlabels = F)

contour(mu, lambda, dens.c, levels = contour4, add = T, drawlabels = F, col = 'blue')
legend("topleft",'(c)',bty = 'n')


#(d) Contours of the optimal factorized approximation--red 
fd <- Vectorize(function(mu, lambda){
  f.d <- sqrt(kappa.optimal / 2/ pi) * exp(-kappa.optimal / 2*(mu - mun)^2) * lambda^(an - 1) * exp(-b.optimal * lambda)
  #logfd=-kappa.optimal*0.5*(mu-mun)^2+(an-1)*log(lambda)-b.optimal*lambda
  
  #return(logfd)
  return(f.d)
})

#logdens.d=outer(mu,lambda,fd)
#dens.d= exp(logdens.d - max(logdens.d))
dens.d <- outer(mu, lambda, fd)

contour5 <- c(5e-43, 1e-40, 1e-38, 1e-37, 1e-35)
contour(mu, lambda, dens, levels = contour1, xlab = 'mu', ylab = 'lambda',col = 'green', drawlabels = F)

contour(mu, lambda, dens.d, levels = contour5, add = T, drawlabels = F, col = 'red')
legend("topleft",'(d)', bty = 'n')

