K <- 2
data <- spirals
data <- as.matrix(data) 

# Kmeans
source('KmeansDemoFaithful.R')	
centers <- data[sample(nrow(data), K),] 
res <- Kmeans(data, centers, euclid, 10)
plot(data[, 1], data[, 2], 
     xlab = "", ylab = "", col = res$clusters)
points(res$centers, col = 1:2, pch = 8, cex = 2)

# Spectral
SpectralCluster <- function(S, k) 
{
  S <- as.matrix(S)
  S <- abs(S)
  D <- diag(rowSums(S))
  V <- eigen(solve(D, S))$vectors[, 2:(k + 1)]
  km <- kmeans(V, k)
  return(km$cluster)
}

S <-  tcrossprod(data)
res  <-  SpectralCluster(S, K)
plot(data[, 1], data[, 2], 
     xlab = "", ylab = "", col = res)