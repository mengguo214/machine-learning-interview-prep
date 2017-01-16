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
# Reference: 
#   http://statistical-research.com/dirichlet-process-infinite-mixture-models-and-clustering/
#   https://github.com/tbroderick/bnp_tutorial/tree/2016mlss

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
    
    # plot cluster assignments in order of appearance
    K = unique(table)
    plot(seq(1, i + 1), table,
         xlab="Sample index", ylab="Cluster by order of appearance",
         xlim=c(0,max(10,i)), ylim=c(0,max(10,length(table))),
         pch=19, main=bquote(rho~"~Dirichlet"  # ~"("~.(a)~",...,"~.(a)~")"
                             ~", K="~.(K))
    )
    
    line <- readline()
    if(line == "x") return("done")
  }
  
}

crp(100, 4)
