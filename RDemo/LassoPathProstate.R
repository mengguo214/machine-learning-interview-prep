# Description: 
#   This script shows Lasso algorithm by using 'prostateStnd' data.
#
# License: GNU v3
#
# Date: October 19, 2016
#
# Authors:
#    Rcode was written by Jiali Lin.
#    Depts. of Statistics, Virginia Tech,
#    Hutcheson Hall, 403K, Blacksburg, VA 24061 

rm(list = ls())

library(R.matlab)
library(glmnet)
library(reshape)
library(ggplot2)
library(latex2exp)

# load data from "pmtk3."
data <- readMat('../data/prostateStnd.mat')
X <- data$Xtrain
Y <- data$ytrain

# fit lasso
model <- glmnet(X, Y)

# extract lambda and estimated coefficient
lambda <- model$lambda
beta.hat <- as.matrix(model$beta)
data <- data.frame(lambda, t(beta.hat))
names(data) <- c("lambda", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")
data <- melt(data, c("lambda"))

# regularization path
ggplot(data, aes(x=lambda,y=value)) +
  geom_point(aes(colour = variable), size = .5) +
  geom_line(aes(colour = variable)) + 
  scale_x_reverse() + 
  coord_cartesian(xlim = c(0, max(lambda) + 2*sd(lambda))) +
  coord_cartesian(ylim = c(min(beta.hat) - 2*sd(beta.hat), 
                                   max(beta.hat) + 2*sd(beta.hat))) +
  xlab(TeX("$\\lambda$")) + 
  ylab(TeX("$\\hat{w}_j(\\lambda)$"))

