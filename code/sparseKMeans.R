#Holding w fixed, update C
UpdateC <- function(x, K, w, C=NULL){
  #num of row
  n <- nrow(x)
  #remove zero-weight features
  x <- x[,w!=0];
  xr <- sweep(x, 2, sqrt(w[w!=0]), "*");
  if(!is.null(C)){
    #num of colomn
    p <- ncol(x)
    #attach cluster information
    xr <- data.frame(data=xr, clu=C)
    #calculate cluster mean
    mu <- aggregate(xr,by=list(xr$clu),FUN=mean)
    #cal new dist including the cluster mean
    xr <- as.matrix(xr[,1:p])
    mu <- mu[,1:p+1]
    mu <- as.matrix(mu)
    distmat <- as.matrix(dist(rbind(xr, mu)))[1:n, (n+1):(n+K)]
    #get the cluster 
    c_n <- apply(distmat, 1, which.min)
    #do the clustering
    if(length(unique(c_n))==K){
      kc <- kmeans(xr,centers=mu)
    }else{
      kc <- kmeans(xr,centers=K, nstart=10)
    }
  }else{
    kc <- kmeans(xr, centers=K, nstart=10)
  }
  return(kc$cluster)
}

FindDelta <- function(a,l1bound){
  #Find optimal delta which statisfy w constrains
  
  #if already statisfies, return 0
  if(l2n(a)==0 || sum(abs(a/l2n(a)))<=l1bound) return(0)
  #binary search lo is 0, hi is maximum absolute element of a
  delta.lo <- 0
  delta.hi <- max(abs(a))-1e-5
  iter <- 1
  while(iter<=15 && (delta.hi-delta.lo)>(1e-4)){
    #cal S(a,delta)
    su <- soft(a,(delta.lo+delta.hi)/2)
   
    if(sum(abs(su/l2n(su)))<l1bound){
      #if stasifies the l1 bound, try a smaller delta
      delta.hi <- (delta.lo+delta.hi)/2
    } else {
      #if not, try last delta
      delta.lo <- (delta.lo+delta.hi)/2
    }
    iter <- iter+1
  }
  return((delta.lo+delta.hi)/2)
}

CalWCSS <- function(x, C, w=NULL){
  wcss.feature <- numeric(ncol(x))
  for(k in unique(C)){
    whichers <- (C==k)
    if(sum(whichers)>1) wcss.feature <- wcss.feature + apply(scale(x[whichers,],center=TRUE, scale=FALSE)^2, 2, sum)
  }
  bcss.feature <- apply(scale(x, center=TRUE, scale=FALSE)^2, 2, sum)-wcss.feature
  if(!is.null(w)) return(list(wcss.feature=wcss.feature, wcss=sum(wcss.feature), wcss.w=sum(wcss.feature*w),
                               bcss.feature=bcss.feature))
  if(is.null(w)) return(list(wcss.feature=wcss.feature, wcss=sum(wcss.feature), bcss.feature=bcss.feature))
}

UpdateW <- function(x, C, l1bound){
  #calculate WCSS
  wcss.feature <- CalWCSS(x, C)$wcss.feature
  #calculate TSS
  tss.feature <- CalWCSS(x, rep(1, nrow(x)))$wcss.feature
  #find delta
  delta <- FindDelta(-wcss.feature+tss.feature, l1bound)
  #soft operation
  wu.unscaled <- soft(-wcss.feature+tss.feature,delta)
  return(wu.unscaled/l2n(wu.unscaled))
}

l2n <- function(vec){
  #calculate l2 norm
  return(sqrt(sum(vec^2)))
}

soft <- function(x,d){
  #soft operation
  return(sign(x)*pmax(0, abs(x)-d))
}


SparseKMeans <- function(x, K, wbounds=NULL, nstart=20){
  #maximum iteration
  maxiter <- 6
  #tolerance
  tol <- 1e-4
  #num of features
  p <- ncol(x);
  #num of instances
  n <- nrow(x)
  #if wbounds is null, initilize it! default 2
  if(is.null(wbounds)) wbounds <- c(3)
  
  C <- kmeans(x, centers=K, nstart=nstart)$cluster
  #Initiliza w as 1/sqrt p
  w <- rep(1/sqrt(p),p) 
  w.old <- rnorm(p)
  
  niter <- 0
  history <- NULL
  while((sum(abs(w-w.old))/sum(abs(w.old)))>1e-4 && niter<maxiter){
    niter <- niter + 1
    w.old <- w
    #update C
    if(niter > 1){
      C <- UpdateC(x, K, w, C);
    }
    #update W
    w <- UpdateW(x,C,wbounds)
    nonzerow <- length(which(w > 1e-10))
    history <- c(history, list(w=w,bcss=sum(CalWCSS(x, C)$bcss.feature*w), nonzerow=nonzerow))
  }
  result <- list(w=w,C=C,nonzerow=length(which(w>1e-10)),history=history)
}

SparseKMeansTuning <- function(x, K, wbounds, nperms=25){
  #num of instances
  n <- nrow(x)
  #num of features
  p <- ncol(x)
  #generate permutation of x
  x.perm <- list()
  for(i in 1:nperms){
    x.perm[[i]] <- matrix(NA,nrow=n, ncol=p)
    for(j in 1:p){
      x.perm[[i]][,j] <- sample(x[,j])
    }
  }
  #initial total sum of s
  totals <- NULL
  #num of nonzero features
  nnonzerows <- NULL
  for(i in 1:length(wbounds)){
    #calculate the totals for each wbound given the original dataset
    km <- SparseKMeans(x, K, wbounds=wbounds[i]);
    nnonzerows <- c(nnonzerows, sum(km$w!=0))
    #calculate BCSS
    bcss <- CalWCSS(x,km$C)$bcss.feature
    totals <- c(totals, sum(km$w*bcss))
  }
  
  permtots <- matrix(NA, nrow=length(wbounds), ncol=nperms)
  for(i in 1:nperms){
    cat('Permuation: ', i, "\n")
    #for each permutation
    for(j in 1:length(wbounds)){
      #for each wbound
      #run sparse K-Means
      perm.km <- SparseKMeans(x.perm[[i]], K, wbounds=wbounds[j]);
      #calculate bcss
      perm.bcss <- CalWCSS(x.perm[[i]],perm.km$C)$bcss.feature
      #calculate weighted sum of bcss
      permtots[j,i] <- sum(perm.km$w*perm.bcss)
    }
  }
  #calculate gaps statistics
  gaps <- (log(totals)-apply(log(permtots),1,mean))
  #calculate gaps standard deviation
  gaps.sd <- apply(log(permtots),1,sd)
  mean.sd <- min(gaps.sd)
  best.index <- which(gaps > max(gaps)-mean.sd)[1]
  #construct the results
  result <- list(gaps=gaps, gaps.sd=gaps.sd,nnonzerows=nnonzerows, best.wbounds=wbounds[best.index], best.windex = best.index)
  return(result)
}