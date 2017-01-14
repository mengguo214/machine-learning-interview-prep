# Description: 
#   This script shows Losses for Binary Classification.
#
# Details:
#    Visualization of 0-1 loss, log loss, hinge loss and exponential loss.
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
#    pmtk3.googlecode.com.

rm(list = ls())

library(ggplot2)
library(latex2exp)

z <- seq(-2, 2, 0.01)
L01 <- as.numeric((sign(z) < 0))
Lhinge <- pmax(0,1 - z)
Lnll <- log2(1+exp(-z))
Lbinom <- log2(1+exp(-2 * z))
Lexp <- exp(-z)

ggplot() + 
  geom_line(aes(x = z,y = L01,color = "0-1")) + 
  geom_line(aes(x = z,y = Lhinge, color = "hinge")) + 
  geom_line(aes(x = z,y = Lnll, color = "logloss")) + 
  xlab(TeX("$\\eta$")) + 
  ylab("loss") + 
  theme(legend.title = element_blank())

ggplot() + 
  geom_line(aes(x = z,y = L01,color = "0-1")) + 
  geom_line(aes(x = z,y = Lexp, color = "exp")) + 
  geom_line(aes(x = z,y = Lnll, color = "logloss")) + 
  xlab(TeX("$\\eta$"))   + 
  ylab("loss") + 
  theme(legend.title = element_blank())
