
# This file contains functions used in simulations for 
# Figure 7 to 10 in the paper.

# The following two are from the best-subset package.
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



sim.xy <- function(n,nval,p,s=5,rho=0,beta.type=1,snr=1){
  
  # function to generate datasets from linear model
  
  # parameter description
  # n         : size of training set
  # nval      : size of testing set
  # p         : number of originla features
  # s         : number of true signal original features
  # rho       : auto-correlation of original features
  # beta.type : pattern of coefficients in the linear model
  # snr       : signal-to-noise ratio
  
  # generate independent features
  x <- matrix(rnorm(n*p),n,p)
  xval <- matrix(rnorm(nval*p),nval,p)
  
  # introduce correlation if needed
  if (rho != 0) {
    inds = 1:p
    Sigma = rho^abs(outer(inds, inds, "-"))
    obj = svd(Sigma)
    Sigma.half = obj$u %*% (sqrt(diag(obj$d))) %*% t(obj$v)
    x = x %*% Sigma.half
    xval = xval %*% Sigma.half
  }else{
    Sigma = diag(1,p)
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
  
  # set the variance of noise based on signal-to-noise ratio
  vmu <- as.numeric(t(beta)%*%Sigma%*%beta)
  sigma <- sqrt(vmu/snr)
  
  # generate response
  y <- as.numeric(x%*%beta+rnorm(n)*sigma)
  yval <- as.numeric(xval%*%beta+rnorm(nval)*sigma)
  
  enlist(x,y,xval,yval,Sigma,beta,sigma)
}



# The following function is used to aggregate results across
# repeated simulations for each setting.

# obj is a list of outputs, each from a single run of the
# wrapper function in the R file Sim_fig_7_to_10.R
agg_data <- function(obj){
  
  n       <- obj[[1]]$n
  nval    <- obj[[1]]$nval
  p       <- obj[[1]]$p
  s       <- obj[[1]]$s
  sigma   <- obj[[1]]$sigma
  SNR     <- obj[[1]]$SNR
  l_SNR   <- length(SNR)
  maxnode <- obj[[1]]$maxnode
  ntree   <- obj[[1]]$ntree
  l_maxnd <- length(maxnode)
  l_ntree <- length(ntree)
  nrep    <- length(obj)
  
  MSE     <- vector(mode = "list", length = l_SNR)
  names(MSE) <- paste('SNR', SNR)
  
  # k: loop for SNRs
  for (k in 1:l_SNR) {
    MSE[[k]]        <- vector(mode = "list", length = 2)
    names(MSE[[k]]) <- c('Train', 'Test')
    
    # i: loop for Train and Test
    for (i in 1:2) {
      MSE[[k]][[i]] <- vector(mode = "list", length = 3)
      names(MSE[[k]][[i]]) <- c('Raw', 'Ave', 'Std')
      
      for (j in 1:3) {
        MSE[[k]][[i]][[j]] <- vector(mode = "list", length = 2)
        names(MSE[[k]][[i]][[j]]) <- c( 'Forest', 'Bagging')
      }
      
      MSE[[k]][[i]][['Raw']][['Forest']] <- matrix(NA, ncol = nrep, nrow = l_maxnd * l_ntree)
      MSE[[k]][[i]][['Raw']][['Bagging']] <- matrix(NA, ncol = nrep, nrow = l_maxnd * l_ntree)
    }
    
    i <- j <- 0
    
    # b: loop for repetitions
    for (b in 1:nrep) {
      MSE[[k]][['Train']][['Raw']][['Forest']][, b]    <- obj[[b]]$MSE_train[[k]]      %>% as.numeric()
      MSE[[k]][['Test']][['Raw']][['Forest']][, b]     <- obj[[b]]$MSE_test[[k]]       %>% as.numeric()
      MSE[[k]][['Train']][['Raw']][['Bagging']][, b]   <- obj[[b]]$MSE_Bag_train[[k]]  %>% as.numeric()
      MSE[[k]][['Test']][['Raw']][['Bagging']][, b]    <- obj[[b]]$MSE_Bag_test[[k]]   %>% as.numeric()
    }
    
    
    for (i in 1:2) {
      
      # j: loop for the default random forest and gagging
      for (j in 1:2) {
        MSE[[k]][[i]][['Ave']][[j]] <- MSE[[k]][[i]][['Raw']][[j]] %>% rowMeans(na.rm = T) %>% matrix(nrow = l_maxnd)
        MSE[[k]][[i]][['Std']][[j]] <- apply(MSE[[k]][[i]][['Raw']][[j]], 1, sd, na.rm = T)/
          sqrt(rowSums(!is.na(MSE[[k]][[i]][['Raw']][[j]]) )  ) %>% matrix(nrow = l_maxnd)
      }
      
      dimnames(MSE[[k]][[i]][['Ave']][['Forest']]) <- dimnames(MSE[[k]][[i]][['Std']][['Forest']]) <-
        dimnames(MSE[[k]][[i]][['Ave']][['Bagging']]) <- dimnames(MSE[[k]][[i]][['Std']][['Bagging']]) <-
        dimnames(obj[[1]]$MSE_test[[1]])
    }
    # The loop for i ends
    
  }
  # The loop for k ends
  
  result <- enlist(n, nval, p, s, sigma, SNR, maxnode, ntree, l_maxnd, l_ntree, nrep, MSE)
  
  return(result)
}
























