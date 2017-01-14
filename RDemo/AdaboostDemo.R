# Description: 
#   Demo of Adaboost using Sonar dataaset
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
#    https://searchcode.com/codesearch/raw/17266536/
# 
# See Also:
#     http://web.stanford.edu/class/stats202/

rm(list = ls())

train <- read.csv("../data/sonar_train.csv", header = FALSE)
test <- read.csv("../data/sonar_test.csv", header = FALSE)
y <- train[, 61]
x <- train[, 1:60]
y.test <- test[, 61]
x.test <- test[, 1:60]
train.error <- rep(0, 50) 
test.error <- rep(0, 50)
f <- rep(0,130) 
f.test <- rep(0,78)
i <- 1
library(rpart)
while(i <= 50){
  w <- exp(-y*f)
  w <- w/sum(w)
  fit <- rpart(y~., x, w, method = "class")
  g <- -1 + 2 * (predict(fit, x)[, 2] > .5) # make -1 or 1
  g.test <- -1 + 2 * (predict(fit, x.test)[, 2] > .5)
  e <- sum(w * (y * g < 0))
  alpha <- .5 * log((1 - e) / e)
  alpha <- 0.1*alpha  #change made for part b
  f <- f + alpha * g
  f.test <- f.test + alpha*g.test
  train.error[i] <- sum(1 * f * y < 0)/130
  test.error[i] <- sum(1 * f.test * y.test < 0)/78
  i <-i+1
}

plot(seq(1, 50), test.error, type = "l", ylim = c(0, .5),
     ylab = "Error Rate",xlab = "Iterations",lwd = 2, main = 'Atul Kumar')
lines(train.error, lwd = 2, col = "purple")
legend(4, .5, c("Training Error", "Test Error"), col = c("purple", "black"), lwd = 2)
