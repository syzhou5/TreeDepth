
# This file contains example codes of simulations for 
# Figure 8 and Table 1.
# Take the data cpu as an example. Simulations on other
# datasets can be carried out similarly.




library(randomForest)
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
source("Functions_fig_8.R")
load("RealData.RData")


# The following is the wrapper function that can be run in
# parallel for repeating the simulation 100 times for each data.
# In each simulation, grow random forests with candidates of mtry 
# and maxnodes, and obtain the cross validation error.

# dat: data to be used
# K  : number of folds for CV

wrapper <- function(nsim, dat, K = 10){
  set.seed(nsim)

  n <- dim(dat)[1]
  p <- dim(dat)[2] -1
  
  
  if (p>10) {
    mtry <- ceiling(seq(10)/10*p)
  }else{
    mtry <- 1:p
  }
  
  if (n>25) { 
    maxnodes <- ceiling(seq(25)/25* floor((1-1/K)*n) )
  }else{
    maxnodes <- 2:floor((1-1/K)*n)
  }
  
  cv.sigma <- cv.rf_noise(dat = dat,  mtry = mtry, maxnodes = maxnodes, K = K)
  
  if(nsim==1){
    result <- enlist(n,p,K,dat,nsim, mtry, maxnodes,cv.sigma)
  }else{
    result <- enlist(cv.sigma)
  }
  
  return(result)
  
}

# Use the wrapper function above to run the simulation 100 times
MSE_out <- lapply(1:100, wrapper, dat = cpu, K=10)

# Aggregate results in the 100 repetitions, draw the boxplot
# of tunability and obtain the optimal value of (mtry, maxnodes)
out <- summary_and_boxplot_real(MSE_out, dt_name = "cpu",
                  mtry0 = 0.3, maxnodes0 = 1)

# boxplot of tunability d 
gp_d <- out$gp_d
# optimal values of (mtry, maxnodes)
out$dat_opt























