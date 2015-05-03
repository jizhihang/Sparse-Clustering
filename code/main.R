source('./sparseKMeans.R')
source('./sparseHierarchical.R')

library(MASS)
set.seed(11)

generateData <- function(p, q, num, bmu){
  K <- 3
  x <- matrix(NA, nrow=num*K, ncol=p)
  
  mu <- c(0,bmu,2*bmu)
  for(i in 1:K){
    bs <- (i-1)*num + 1
    be <- (i-1)*num + num
    for(j in 1 : p){
      if(j <= q){
        x[bs:be,j] = rnorm(num, mean=mu[i])
      }else{
        x[bs:be,j] = rnorm(num)
      }
    }
  }
  return(x)
}
set.seed(111)
x <- generateData(100,50,20,0.6)
#Sparse K-means
k.perm <- SparseKMeansTuning(x,K=3,wbounds=seq(3,10,len=15), nperm=10)
best.index = k.perm$best.windex
#plot gap statistics
plot(k.perm$gaps,main="Gap Statistics", type='b',xlab="s",ylab="gap",col = "dark blue", lwd = 4)
points(best.index,k.perm$gaps[best.index],col="red",type = 'p',pch = 16)
#plot non-zero features
plot(k.perm$nnonzerows,main="Non-Zero Features",type="b",xlab="s",ylab="num of nonzero feature", col = "dark blue", lwd = 4)
points(best.index,k.perm$nnonzerows[best.index],col="red",type = 'p',pch = 16)

par(mfrow = c(2,2),pty='s')
k.out <- SparseKMeans(x,K=3,wbounds=3)
plot(k.out$w, col="dark blue", pch=16, type="p", xlab="Feature Index", ylab="weight", main="s=3")

k.out <- SparseKMeans(x,K=3,wbounds=5)
plot(k.out$w, col="dark blue", pch=16, type="p", xlab="Feature Index", ylab="weight", main="s=5")

k.out <- SparseKMeans(x,K=3,wbounds=10)
plot(k.out$w, col="dark blue", pch=16, type="p", xlab="Feature Index", ylab="weight", main="s=10")

k.out <- SparseKMeans(x,K=3,wbounds=20)
plot(k.out$w, col="dark blue", pch=16, type="p", xlab="Feature Index", ylab="weight", main="s=100")

#w
pa <- seq(100,1000,100)
mu <- 0.6

rndk <- rep(0,length(pa))
rndkm <- rep(0,length(pa))
for(i in 1 : length(pa)){
  ks <- 0
  kms <- 0
  print(i)
  for(j in 1 : 20){
    x <- generateData(pa[i],50,20,mu)
    kc <- SparseKMeans(x,K=3,wbounds=5)$C
    kmc <- kmeans(x,centers=K, nstart=20)$cluster
    true.id <- c(rep(1,20),rep(2,20),rep(3,20))
   
    ks <- ks + as.numeric(RRand(true.id,kc)[1])
    kms <- kms + as.numeric(RRand(true.id,kmc)[1])
  }
  rndk[i] <- ks/20
  rndkm[i] <- kms/20
}
plot(pa,rndkm, type='b',xlab="p",ylab="Rand Index",col = "dark blue", lwd = 4, main="Rand Index")
lines(pa,rndk, type='b',xlab="p",ylab="Rand Index",col = "red", lwd = 4)
legend("topright", c("Sparse K-Means","Standard K-Means"), col = c("red","dark blue"),lty=1,lwd=4)

h.out <- sparseHierarchicalTuning(x,wbound=c(1.5,2:6))
