
# This file contains example codes of simulations for 
# Figure 7 to 10.

# The following is for Figure 7. To generate the other 
# 3 plots, just modify the parameters for data generation.



# In the following, several possible values of ntree are
# under consideration, while only results corresponding 
# to the largest value of ntree, which is 500, are plotted.


library(randomForest)
library(tidyverse)
require(ggplot2)
require(grid)
require(gridExtra)


# The following is the wrapper function that can be run in
# parallel for repeating the simulation 100 times.


wrapper <- function(nsim){
  
  set.seed(nsim)
  
  # Specify parameters for generating data
  n            <- 100
  nval         <- 100
  p            <- 10
  s            <- 5
  rho          <- 0.35
  beta.type    <- 2
  SNR          <- c(0.05, 0.09, 0.14, 0.25, 0.42, 0.71, 
                    1.22, 2.07, 3.52, 6)
  K            <- length(SNR)
  
  # Parameters for growing forests
  maxnode      <- c(seq(2, 8, by = 2),
                    seq(10, n/2, by = 5),
                    seq(n/2 + 25, n, by = 25))
  l_maxnode    <- length(maxnode)
  ntree        <- c(1, seq(10, 100, by = 10), seq(150, 500, by = 50))
  l_ntree      <- length(ntree)
  
  # Initialize lists/vectors to store results
  MSE_train    <- MSE_test    <- vector(mode = "list", length = K)
  MSE_Bag_train  <- MSE_Bag_test  <- vector(mode = "list", length = K)
  sigma        <- numeric(K)
  
  for (k in 1:K) {
    cat(k,"\n")
    MSE_train[[k]]  <- MSE_test[[k]] <- 
      MSE_Bag_train[[k]] <- MSE_Bag_test[[k]] <- 
      matrix(NA, nrow = l_maxnode, ncol = l_ntree)
    
    rownames(MSE_train[[k]]) <- rownames(MSE_test[[k]]) <- 
      rownames(MSE_Bag_train[[k]]) <- rownames(MSE_Bag_test[[k]]) <- 
      paste("maxnode", maxnode)
    colnames(MSE_train[[k]]) <- colnames(MSE_test[[k]]) <- 
      colnames(MSE_Bag_train[[k]]) <- colnames(MSE_Bag_test[[k]]) <- 
      paste("ntree", ntree)
    
    # Generate Data
    xy.obj     <- sim.xy(n, nval, p, s, rho, beta.type, SNR[k])
    sigma[k]   <- xy.obj[['sigma']]
    
    x          <- xy.obj[["x"]]
    y          <- xy.obj[["y"]] 
    xval       <- xy.obj[["xval"]]
    yval       <- xy.obj[["yval"]]
    
    for (i in 1:l_maxnode) {
      
      for (j in 1:l_ntree) {
        rf <- randomForest(x = x, y = y, xtest = xval, ytest = yval,
                           ntree = ntree[j], maxnodes = maxnode[i],
                           nodesize = 1)
        MSE_train[[k]][i,j] <- rf[["mse"]][ntree[j]]
        MSE_test[[k]][i,j]  <- rf[["test"]][["mse"]][ntree[j]]
        
        bag <- randomForest(x = x, y = y, xtest = xval, ytest = yval,
                            ntree = ntree[j], maxnodes = maxnode[i],
                            nodesize = 1, mtry = p)
        MSE_Bag_train[[k]][i,j] <- bag[["mse"]][ntree[j]]
        MSE_Bag_test[[k]][i,j]  <- bag[["test"]][["mse"]][ntree[j]]
      }
    }
  }
  
  names(MSE_train) <- names(MSE_test) <- 
    names(MSE_Bag_train)  <- names(MSE_Bag_test)  <-
    paste("SNR", SNR)
  
  if (nsim==1) {
    result <- enlist(n, nval, p, s, rho, beta.type, SNR, sigma, 
                     l_maxnode, l_ntree, maxnode, ntree,
                     MSE_train, MSE_test,
                     MSE_Bag_train, MSE_Bag_test)
  }else{
    result <- enlist(MSE_train, MSE_test,
                     MSE_Bag_train, MSE_Bag_test)
  } 
  return(result)
}


# Use the wrapper function above to run the simulation 100 times
MSE_out <- lapply(1:100, wrapper)

# Aggregate results in the 100 repetitions

dt_summary <- agg_data(obj = MSE_out)

maxnodes <- dt_summary$maxnode
ntree    <- dt_summary$ntree
SNR      <- dt_summary$SNR
l_ntree  <- length(ntree)

# the number of SNRs under consideration
ls       <- length(dt_summary$SNR)
gp_rfbag <- vector(mode = "list", length = ls)

for (j in 1:ls) {
  
  dat_rfbag <- data.frame(maxnodes = rep(maxnodes, 2),
                          MSE      = c(
                            dt_summary$MSE[[j]]$Test$Ave$Forest[, l_ntree] %>% as.numeric(), 
                            dt_summary$MSE[[j]]$Test$Ave$Bagging[, l_ntree] %>% as.numeric()
                          ),
                          sd       = c(
                            dt_summary$MSE[[j]]$Test$Std$Forest[, l_ntree] %>% as.numeric(), 
                            dt_summary$MSE[[j]]$Test$Std$Bagging[, l_ntree] %>% as.numeric()
                          ),
                          Model    = rep(c( "Random Forest", "Bagging"), 
                                         each = length(maxnodes)
                          ) %>% factor()
                          )
  
  gp_rfbag[[j]]  <- dat_rfbag %>% ggplot( aes(x = maxnodes, y = MSE, color = Model) ) +
    geom_point(aes(shape = Model) ) +
    geom_line(lwd = 1) +
    geom_errorbar(aes(ymin = MSE - sd, ymax = MSE + sd)) +
    scale_color_manual(values = c("lightskyblue", "red")) +
    xlab("Maxnodes") +
    ylab("Test MSE") + 
    ggtitle(paste0("SNR = ", SNR[j])) +
    theme(legend.position = c(0.98, 0.98),
          legend.justification = c(1,1),
          plot.title = element_text(hjust = 0.5,size = rel(2.5), face = 'bold'),
          legend.key.size = unit(0.04, "npc"),
          legend.key = element_rect(fill = "gray88"),
          legend.spacing.y = unit(.01, "npc"),
          legend.text = element_text(size = rel(2)),
          legend.title = element_text(size = rel(2)),
          axis.title = element_text(hjust = 0.5,size = rel(2)),
          axis.text  = element_text(size = rel(1.5))
    )
  
  if (j!= ls) {
    gp_rfbag[[j]] <- gp_rfbag[[j]] + theme(legend.position = "None")
  }
}

gp_rfbag_out <- grid.arrange(grobs = gp_rfbag, ncol = 2)
gp_rfbag_out






















