# Description: 
#   Illustration of the k-Means algorithm using the Old Faithful dataset.
#
# Details:
#    K: the number of clusters is set to 2.
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
#    http://stackoverflow.com/questions/31571236/my-own-k-means-algorithm-in-r

rm(list = ls())
 
# The K means algorithm that uses the Euclidean distance
euclid <- function(points1, points2) {
  distanceMatrix <- matrix(NA, nrow = dim(points1)[1], ncol = dim(points2)[1])
  for(i in 1:nrow(points2)) {
    distanceMatrix[, i] <- sqrt(rowSums(t(t(points1)-points2[i, ])^2))
  }
  distanceMatrix
}

# K-Means algorithm
Kmeans <- function(x, centers, dist.fun, nItter) {
  cluster.history <- vector(nItter, mode = "list")
  center.history <- vector(nItter, mode = "list")
  
  for(i in 1:nItter) {
    distsToCenters <- dist.fun(x, centers)
    clusters <- apply(distsToCenters, 1, which.min)
    centers <- apply(x, 2, tapply, clusters, mean)
    
    # Saving history
    cluster.history[[i]] <- clusters
    center.history[[i]] <- centers
  }
  list(clusters = cluster.history[[nItter]], centers = center.history[[nItter]])
}

# Demo using Old Faithful dataset
K <- 2
data <- faithful
data <- as.matrix(data) 
centers <- data[sample(nrow(data), K),] 
res <- Kmeans(data, centers, euclid, 10)
plot(faithful$eruptions, faithful$waiting, 
     xlab = "eruptions", ylab = "waiting",
     col = res$clusters)
points(res$centers, col = 1:2, pch = 8, cex = 2)


