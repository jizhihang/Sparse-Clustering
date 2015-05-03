library(sparcl)
library(MASS)
set.seed(111)

num <- 300
mu1 <- c(0,5)
sigma <- matrix(c(20,5,5,2),2,2)
x1 <- mvrnorm(num, mu1, sigma)
mu2 <- c(20,5)
x2 <- mvrnorm(num, mu2, sigma)
x <- rbind(x1,x2)
x <- scale(x, TRUE, TRUE)

par(mfrow = c(2,1),pty='s')

#K-means
cl <- kmeans(x,2)
plot(x,col=cl$cluster, main='Standard K-means')

#Sparse K-means
km.perm <- KMeansSparseCluster.permute(x,K=2,wbounds=seq(3,7,len=15),nperms=5)
km.out <- KMeansSparseCluster(x,K=2,wbounds=km.perm$bestw)
plot(x,col=km.out[[1]]$Cs, main='Sparse K-means')
