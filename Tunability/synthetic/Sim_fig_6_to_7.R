
# This file contains example codes of simulations for 
# Figure 6, 7 and B7.

# The following is for the Low setting. For the other 
# settings, just modify the parameters for data generation.





library(randomForest)
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
source("Functions_fig_6_to_7.R")


# The following is the wrapper function that can be run in
# parallel for repeating the simulation 100 times.
# In each simulation, grow random forests with the candidates
# of mtry and maxnodes and obtain their training and test MSEs.


wrapper <- function(nsim){
  
  set.seed(nsim)
  
  # Specify parameters for generating data
  n            <- 100
  nval         <- 100
  p            <- 10
  s            <- 5
  rho          <- 0.35
  beta.type    <- 2
  SNR          <- c(0.05, 0.09, 0.14, 0.25, 0.42, 0.71, 1.22, 2.07, 3.52, 6)
  
  K            <- length(SNR)
  
  # Specify parameters for tree construction
  maxnode      <- c(seq(2, 8, by = 2),
                    seq(10, n/2, by = 5),
                    seq(n/2 + 25, n, by = 25))
  l_maxnode    <- length(maxnode)
  
  mtry         <- seq(10)/10
  l_mtry       <- length(mtry)
  
  ntree        <- 500
  
  
  # Initialize lists/vectors to store results
  MSE_train    <- MSE_test    <- vector(mode = "list", length = K)
  sigma        <- numeric(K)
  
  for (k in 1:K) {
    MSE_train[[k]]  <- MSE_test[[k]] <-
      matrix(NA, nrow = l_maxnode, ncol = l_mtry)
    
    rownames(MSE_train[[k]]) <- rownames(MSE_test[[k]]) <-
      paste("maxnode", maxnode)
    colnames(MSE_train[[k]]) <- colnames(MSE_test[[k]]) <-
      paste("mtry", mtry)
    
    xy.obj     <- sim.xy(n, nval, p, s, rho, beta.type, SNR[k])
    sigma[k]   <- xy.obj[['sigma']]
    
    x          <- xy.obj[["x"]]
    y          <- xy.obj[["y"]] 
    xval       <- xy.obj[["xval"]]
    yval       <- xy.obj[["yval"]]
    
    
    for (i in 1:l_maxnode) {
      # loop for tree depth
      
      for (j in 1:l_mtry) {
        # loop for randomness
        cat(k, i, j, "\n")
        rf <- randomForest(x = x, y = y, xtest = xval, ytest = yval,
                           ntree    = ntree, 
                           maxnodes = maxnode[i],
                           nodesize = 1,
                           mtry     = ceiling( mtry[j]*p))
        
        MSE_train[[k]][i,j] <- rf[["mse"]][ntree]
        MSE_test[[k]][i,j]  <- rf[["test"]][["mse"]][ntree]
        
      }# end for randomness
      
    }# end for depth
    
    
  }# end for SNR
  
  names(MSE_train) <- names(MSE_test) <- 
    paste("SNR", SNR)
  
  
  if (nsim==1) {
    result <- enlist(n, nval, p, s, rho, beta.type, SNR, sigma, 
                     l_maxnode, maxnode, ntree, mtry, l_mtry,
                     MSE_train, MSE_test)
  }else{
    result <- enlist(MSE_train, MSE_test)
  } 
  return(result)
}


# Use the wrapper function above to run the simulation 100 times
MSE_out <- lapply(1:100, wrapper)

# Aggregate results in the 100 repetitions and draw the boxplot
# and scatterplot of tunability and the plot of optimal values
# of (mtry, maxnodes)
out <- summary_and_plot(MSE_out, 
                 g_title = "Low", mtry0 = 0.3, maxnode0 = 1)
# boxplot of tunability d 
gp_d <- out$gp_d
# scatterplot of tunability d
gp_d_ave <- out$gp_d_ave
# optimal values of (mtry, maxnodes)
gp_opt <- out$gp_opt























