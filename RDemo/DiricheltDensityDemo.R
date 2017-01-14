# Description: 
#   This script shows Density plot of Dirichlet distribution
#
# Details:
#    Density plot of Dirichlet(2,2,2)
#
# License: GNU v3
#
# Date: October 19, 2016
#
# Authors:
#    Rcode was modified by Jiali Lin.
#    Depts. of Statistics, Virginia Tech,
#    Hutcheson Hall, 403K, Blacksburg, VA 24061 
#
# Source: 
#   http://www.apps.stat.vt.edu/zhu/teaching/2016/6474/6474_2016.htm

rm(list = ls())

library(gtools)
alpha = c(2,2,2)
p1 = seq(0,1,length=50)
p2 = p1
dens = matrix(0, 50, 50)

for (i in 1:50)
  for (j in 1:50)
    dens[i,j] = ifelse(p1[i] + p2[j] > 1, NA, 
                       ddirichlet(c(p1[i], p2[j], 1 - p1[i] - p2[j]), alpha = c(2, 2, 2)))

image(x = p1,y = p2, z = dens, 
      xlim = c(0,1),ylim = c(0,1), 
      zlim = c(min(dens,na.rm = T), max(dens, na.rm = T )),
      col = heat.colors(30),
      main = "Contour plot of the density of Dirichlet(2,2,2)",
      xlab = "p1", ylab = "p2")
contour(x = p1,y = p2,z = dens, add = TRUE, drawlabels = FALSE, col = "grey", lty = 2)
