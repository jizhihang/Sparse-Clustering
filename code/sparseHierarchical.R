CalUW <- function(distmat, wbound, maxiter=15){
  #num of features
  p = ncol(distmat)
  #intialize u and w
  u <- rnorm(nrow(distmat))
  w <- rep(1, p)
  w.old <- rnorm(p)
  niter <- 1
  while(niter <= maxiter && (sum(abs(w.old-w))/sum(abs(w.old)))>1e-4){
      niter <- niter + 1
      w.old <- w
      #u = Dw/||Dw||^2
      u <- distmat%*%matrix(w,ncol=1)
      u <- u/l2n(u)
      #a = Du
      a <- matrix(u,nrow=1)%*%distmat
      #cal delta
      delta <- FindDelta(a,wbound)
      #w = S(a+,delta)/||S(a+,delta)||^2
      w <- soft(a,delta)
      w <- w/l2n(w)
  }
  u <- distmat%*%matrix(w,ncol=1)/sum(w)
  u <- u/l2n(u)
  w <- w/l2n(w)
  #calculate object value u*D*w
  crit <- sum(u*(distmat%*%matrix(w,ncol=1)))
  #restore distance matrix
  u2 <- matrix(0,nrow=ceiling(sqrt(2*length(u))),ncol=ceiling(sqrt(2*length(u))))
  u2[lower.tri(u2)] <- u
  u <- as.matrix(as.dist(u2))/sqrt(2)
  return(list(u=u, w=w, crit=crit))
}

sparseHierarchical <- function(x,wbound=NULL,method="complete", maxiter=15){
  #Only support squared distance
  
  #calculate distance matrix
  distmat <- matrix(distfun(x), ncol=ncol(x))
  distmat <- distmat^2
  #initial wbound
  if(is.null(wbound)) wbound <- .5*sqrt(ncol(distmat))
  #calculate u and w
  out <- CalUW(distmat, wbound, maxiter = maxiter)
  out <- list(hc = hclust(as.dist(out$u), method = method), w = out$w, u = out$u, 
              crit = out$crit, dists = distmat, wbound = wbound)
  return(out)
}

sparseHierarchicalTuning <- function(x, wbounds, nperms=10){
  #initialize
  tots <- rep(NA,length(wbounds))
  permtots <- matrix(NA, nrow=length(wbounds), ncol=nperms)
  nnonzerows <- rep(NA,length(wbounds))
  #for each wbound, run sparse hierarchical clustering on original data set
  for(i in 1:length(wbounds)){ 
    out <- sparseHierarchical(x,  wbound=wbounds[i])
    #store the results
    nnonzerows[i] <- sum(out$w!=0)
    tots[i] <- out$crit
  }
  
  permdists <- out$dists
  #permutation n times
  for(k in 1:nperms){
    for(j in 1:ncol(permdists)) permdists[,j] <- sample(permdists[,j])
    for(i in 1:length(wbounds)){
      #run hierarchical clustering
      perm.out <- sparseHierarchical(x,wbound=wbounds[i])
      permtots[i,k] <- max(perm.out$crit)
    }
  }
  gaps <- (log(tots)-apply(log(permtots),1,mean))
  out <- list(tots=tots, permtots=permtots, nnonzerows=nnonzerows, gaps=gaps, sdgaps=apply(log(permtots),1,sd), wbounds=wbounds, best.wbound=wbounds[which.max(gaps)], dists=out$dists)
  return(out)
}