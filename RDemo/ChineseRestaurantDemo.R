# Description: 
#   This script shows an example of chinese restaurant process.
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
#   http://statistical-research.com/dirichlet-process-infinite-mixture-models-and-clustering/

rm(list = ls())

crp <- function(num.customers, alpha) {
  table <- 1
  next.table <- 2
  for (i in 1:num.customers) {
    if (runif(1, 0, 1) < alpha / (alpha + i)) {
      # Add a new ball color.
      table <- c(table, next.table)
      next.table <- next.table + 1
    } else {
      # Pick out a ball from the urn, and add back a
      # ball of the same color.
      select.table <- table[sample(1:length(table), 1)]
      table <- c(table, select.table)
    }
  }
  table
}

crp(100, 4)

plot(table(crp(100, 2)), 
     xlab="Table Index", 
     ylab="Frequency")
