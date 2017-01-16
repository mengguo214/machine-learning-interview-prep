# This script shows 
#   an example of polya urn process.
#
# Desciption:
#   Assume F ~ DP(alpha, G0) with alpha = 5, G0 ~ N (0, 1). 
#   Let theta.i ~ F, generate 10000 samples from F using the Polya urn scheme.
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
# Source:
#   http://www.apps.stat.vt.edu/zhu/teaching/2016/6474/6474_2016.htm

polya_urn_model <- function(G_0, N, alpha) {
  balls <- c()
  
  for (i in 1:N) {
    if (runif(1) < alpha / (alpha + length(balls))) {
      # Add a new ball color.
      new_color <- G_0()
      balls <- c(balls, new_color)
      } 
    else {# Pick out a ball from the urn, and add back a
      # ball of the same color.
      ball <- balls[sample(1:length(balls), 1)]
      balls <- c(balls, ball)
      }
  }
  return(balls)
}

G_0 <- function(){
  rnorm(1, mean = 0, sd = 1)
}
N = 100;  alpha = 2;

theta <- polya_urn_model(G_0, N, alpha)

par(mfrow = c(1, 2))
hist(theta, breaks = 20, xlim = c(-1, 1), col = "black")

plot(theta, seq(1, N),
     ylab = "Sample index", xlab = "Cluster by order of appearance",
     ylim = c(0,max(10,N)), xlim = c(-1, 1), pch = 19)
