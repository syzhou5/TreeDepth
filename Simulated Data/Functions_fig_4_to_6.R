
# This file contains functions used in simulations for 
# Figure 4 to 6 in the paper.


# enlist is from the best-subset package.
# https://github.com/ryantibs/best-subset
enlist <- function (...) 
{
  result <- list(...)
  if ((nargs() == 1) & is.character(n <- result[[1]])) {
    result <- as.list(seq(n))
    names(result) <- n
    for (i in n) result[[i]] <- get(i)
  }
  else {
    n <- sys.call()
    n <- as.character(n)[-1]
    if (!is.null(n2 <- names(result))) {
      which <- n2 != ""
      n[which] <- n2[which]
    }
    names(result) <- n
  }
  result
}


# The following sim.xy2 is used to generate data from linear models 
# considered in Hastie et al. 2020.
# It is almost the same as the sim.xy function in the best-subset
# package except that an extra test set is also generated.

sim.xy2 <- function(n,p,nval,ntest,rho=0,s=5,beta.type=1,snr=1){
  # generate covariates
  x <- matrix(rnorm(n*p),n,p)
  xval <- matrix(rnorm(nval*p),nval,p)
  xtest <- matrix(rnorm(ntest*p),ntest,p)
  
  if (rho!=0) {
    ind <- 1:p
    Sigma <- rho^abs(outer(ind,ind,FUN = "-"))
    obj <- svd(Sigma)
    Sigma.half <- obj$u %*% sqrt(diag(obj$d)) %*% t(obj$v)
    x <- x%*% Sigma.half
    xval <- xval%*% Sigma.half
    xtest <- xtest%*% Sigma.half
  }else{
    Sigma <- diag(1,p)
  }
  
  # Generate underlying coefficients
  s = min(s,p)
  beta = rep(0,p)
  if (beta.type==1) {
    beta[round(seq(1,p,length=s))] = 1
  } else if (beta.type==2) {
    beta[1:s] = 1
  } else if (beta.type==3) {
    beta[1:s] = seq(10,0.5,length=s)
  } else if (beta.type==4) {
    beta[1:6] = c(-10,-6,-2,2,6,10)
  } else {
    beta[1:s] = 1
    beta[(s+1):p] = 0.5^(1:(p-s))
  }
  
  # Set snr based on sample variance on infinitely large test set
  vmu <- as.numeric(t(beta)%*%Sigma%*%beta)
  sigma <- sqrt(vmu/snr)
  
  # Generate response
  y <- as.numeric(x%*%beta+rnorm(n)*sigma)
  yval <- as.numeric(xval%*%beta+rnorm(n)*sigma)
  ytest <- as.numeric(xtest%*%beta+rnorm(n)*sigma)
  
  
  enlist(x,y,xval,yval,xtest,ytest,Sigma,beta,sigma,n,p,nval,ntest)
}



# The following is used to contruct random forests and evaluate their
# performance on the validation and test sets.

sim.rf.optTree <- function(n, p, nval, ntest, rho=0, s=5, beta.type=1, snr=1,
                           ndsize, ntree=500, replace=T, 
                           sampsize=ifelse(replace,n,ceiling(0.632*n)),
                           mtry=max(floor(p/3),1)){
  ##### Description of Parameters #####
  ### Parameters for generating data
  # n: sample size of training set
  # p: the number of features
  # nval: sample size of validation set
  # ntest: sample size of testing set
  # rho: the parameter determining the pairwise correlation
  # s: the number of nonzero coefficients
  # beta.type: determining the pattern of beta
  # snr: a vector of signal-to-noise ratio whose length is l.s
  
  ### Parameters for the random forest
  # ndsize: node size to  for each tree
  # ntree: the number of trees to be built
  # replace: bootstrap(T) or subsample(F)
  # sampsize: sample size of bootstrapped sample/subsample
  # mtry: the number of candidates at each split. 
  
  
  nnode <- length(ndsize)
  l.s <- length(snr)
  
  # Initiate lists with matrices as elements for storing values
  # Each row of a matrix corresponds to a tree for a specific snr
  muhat.val.tree <-  muhat.test.tree <- vector(mode="list",length = l.s)
  names(muhat.val.tree) <- names(muhat.test.tree) <- apply(as.matrix(round(snr,2)),2,function(x){paste("snr=",x,sep="")})
  
  muhat.test <- matrix(NA,ntest,l.s)
  mse.test <-  numeric(l.s)
  
  sigma <- numeric(l.s)
  ndsize.tune <- numeric(l.s)
  
  
  # Iterate over snr
  for (i in 1:l.s) {
    
    muhat.test.tree[[i]]<- vector(mode="list",length=nnode)
    muhat.val.tree[[i]] <- vector(mode="list",length=nnode)
    
    for (j in 1:nnode) {
      muhat.test.tree[[i]][[j]] <- matrix(NA,ntest,ntree)
      muhat.val.tree[[i]][[j]] <- matrix(NA,nval,ntree)
    }
    
    # Generate datasets
    xy.obj <- sim.xy2(n,p,nval,ntest,rho,s,beta.type,snr=snr[i])
    sigma[i] <- xy.obj$sigma
    
    # Build trees
    for (b in 1:ntree) {
      
      # print(c(b,i))
      # Draw a bootstrapped sample/subsample such that all the following trees 
      # with different nodesizes use the sample training sample
      train <- sample(1:n,sampsize,replace=replace)
      x.t <- xy.obj$x[train,,drop=F]
      y.t <- xy.obj$y[train,drop=F]
      
      # Set replece=F to make sure the sample used is the bootstrapped sample we draw above 
      # so that rfs with different nodesizes use the same sample
      
      for (j in 1:nnode) {
        rf<- randomForest(x=x.t,y=y.t,xtest=xy.obj$xval,ytest=xy.obj$yval,ntree=1,
                               replace=F,sampsize=sampsize,nodesize=ndsize[j],mtry=mtry,keep.forest = T)
        
        muhat.val.tree[[i]][[j]][,b] <- rf$test$predicted
        muhat.test.tree[[i]][[j]][,b] <- predict(rf,newdata=xy.obj$xtest)
      }
    }
    
    # Aggregate and tune the forest
    # store the prediction and mse of the forest on the validation/test set for each nodesize
    muhat.val.i.temp <- muhat.test.i.temp <- vector(mode="list",length=nnode)
    mse.val.i.temp <- mse.test.i.temp <- numeric(nnode)
    
    for (j in 1:nnode) {
      
      muhat.val.i.temp [[j]] <- as.matrix(rowMeans(muhat.val.tree[[i]][[j]]))
      muhat.test.i.temp [[j]] <- as.matrix(rowMeans(muhat.test.tree[[i]][[j]]))
      
      mse.val.i.temp[j] <- mean((muhat.val.i.temp[[j]] - xy.obj$yval)^2)
      mse.test.i.temp[j] <- mean((muhat.test.i.temp[[j]] - xy.obj$ytest)^2)
    }
    
    ind.tune <- which.min(mse.val.i.temp)
    ndsize.tune[i] <- ndsize[ind.tune]
    muhat.test[,i] <- muhat.test.i.temp[[ind.tune]]
    mse.test[i] <- mse.test.i.temp[ind.tune]
    
  }
  
  enlist(mse.test,muhat.test, ndsize.tune,sigma)
  
}

