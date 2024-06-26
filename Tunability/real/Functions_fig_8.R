
# This file contains functions used in simulations for 
# Figure 8 in the paper.

# The following function enlist is from the best-subset package.
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


se <- function(x, margin =2){
  
  if (margin ==2) {
    out <- apply(x, margin, sd, na.rm = T)/
      sqrt(colSums(!is.na(x)))
  }else if(margin ==1){
    out <- apply(x, margin, sd, na.rm = T)/
      sqrt(rowSums(!is.na(x)))
  }
  
  
  return(out)
}




# K is the number of folds for CV

cv.rf_noise <- function(dat, mtry, maxnodes, K=10, ntree = 500){
  if (K <3) {
    stop("K must be bigger than 2; K=10 recommended")
  }
  
  n <- dim(dat)[1]
  p <- dim(dat)[2]-1
  
  
  l_mtry <- length(mtry)
  l_maxnodes <- length(maxnodes)
  
  foldid  <- sample(rep(seq(K), length = n))
  
  cvm <- matrix(NA, nrow = l_mtry, ncol = l_maxnodes)
  
  for (i in 1:l_mtry) {
    
    for (j in 1:l_maxnodes) {
      
      cv_error <- numeric(K)
      
      for (f in 1:K) {
        which = foldid == f
        
        rf <- randomForest(x = dat[-which, 1:p],
                           y = dat[-which, p+1],
                           xtest = dat[which, 1:p],
                           ytest = dat[which, p+1],
                           mtry = mtry[i],
                           nodesize = 1,
                           maxnodes = maxnodes[j],
                           ntree = ntree
        )
        
        cv_error[f] <- rf$test$mse[ntree]
      } # end for K-fold
      
      cvm[i,j] <- mean(cv_error)
      
    } # end for depth maxnodes
    
  } # end for randomness mtry
  
  rownames(cvm) <- paste0("mtry = ", mtry)
  colnames(cvm) <- paste0("maxnodes = ", maxnodes)
  
  return(cvm)
  
}





# The following function is used to aggregate results across
# repeated simulations for each dataset.

# FinalResult : outputs of simulations from running the wrapper
# dt_name    : name of the dataset; included in the output for 
#              reference and used as the title of plots
# mtry0       : default for mtry, as a proprotion of p
# maxnode0    : default for maxnodes, as a proportion of n
# lgd_pos     : position of legends
# ymax_box    : upper limits for y-axis in boxplots


summary_and_boxplot_real <- function(FinalResult, 
                                     dt_name,
                                     mtry0=0.3, maxnodes0=1,
                                     lgd_pos = c(0.02, 0.98),
                                     ymax_box = NULL,
                                     width = 8,
                                     height = 6,
                                     device = "pdf",
                                     save = T){
  
  n <- FinalResult[[1]]$n
  p <- FinalResult[[1]]$p
  K <- FinalResult[[1]]$K
  dat <- FinalResult[[1]]$dat
  
  mtry <- FinalResult[[1]]$mtry
  maxnodes <- FinalResult[[1]]$maxnodes
  var_y <- FinalResult[[1]]$dat[,(p+1)] %>% var()
  nrep <- length(FinalResult)
  
  l_maxnodes <- length(maxnodes)
  l_mtry <- length(mtry)
  
  # Calculate the normalizing constants for tunability
  normal_cst <- var_y
  
  # Initialize lists for tunability d
  d_RF <- d_mtry <- d_depth <- d_add <- numeric(nrep)
  
  
  opt_RF  <- matrix(NA, nrow = nrep, ncol = 2)
  colnames(opt_RF) <- c("maxnodes", "mtry")
  
  # Index of default values of parameters
  mtry_dft <- which.min(abs(mtry - ceiling(mtry0*p)))
  depth_dft <- which.min(abs(maxnodes - maxnodes0*n))

  for (i in 1:nrep) {
    RTE     <- FinalResult[[i]]$cv.sigma/normal_cst
    RTE_dft <- RTE[mtry_dft, depth_dft]
    
    # optimal combination of mtry and maxnodes
    ind_opt_RF <- which(RTE == min(RTE), arr.ind = T)
    opt_RF[i, ] <- c(maxnodes[ind_opt_RF[2]]/n, mtry[ind_opt_RF[1]]/p)
    
    d_RF[i]    <- RTE_dft - min(RTE)
    d_mtry[i]  <- RTE_dft - min(RTE[, depth_dft])
    d_depth[i] <- RTE_dft - min(RTE[mtry_dft, ])
    d_add[i]   <- min(c(min(RTE[, depth_dft]),
                        min(RTE[mtry_dft, ]))) - min(RTE)
  }    

  
  opt_RF_med  <- apply(opt_RF, 2, median, na.rm=TRUE)%>% matrix(nrow = 2)
  opt_RF_mad   <- (apply(opt_RF, 2, mad, na.rm=TRUE) /
                     sqrt(colSums(!is.na(opt_RF)))) %>% matrix(nrow = 2)
  rownames(opt_RF_med) <- rownames(opt_RF_mad) <- c("maxnodes","mtry")
  
  # Generate boxplots of tunability d
  dat_d <- data.frame(d         = c(d_RF    %>% as.numeric(),
                                    d_mtry  %>% as.numeric(),
                                    d_depth %>% as.numeric(),
                                    d_add   %>% as.numeric()),
                      Type      = c("RF", "mtry", "maxnodes", "Additional") %>% 
                                    rep(each = nrep) %>% as.factor(),
                      mtry_dft  = mtry[mtry_dft] %>% as.factor(),
                      depth_dft = maxnodes[depth_dft] %>% as.factor(),
                      dataset   = dt_name %>% as.factor()
  )
  
  
  gp_d <- ggplot(dat_d, aes(x = Type, y = d)) +
    geom_boxplot(aes(color = Type)) +
    xlab("Type") +
    ylab("Tunability") + 
    ggtitle(dt_name)
  
  if (!is.null(ymax_box)) {
    gp_d <- gp_d + coord_cartesian(ylim = c(0,ymax_box)) 
  }
  
  gp_d <- gp_d + 
    theme_bw() +
    theme(legend.position = "None",
          legend.justification = c(0,1),
          plot.title = element_text(hjust = 0.5,size = rel(2.5), face = 'bold'),
          legend.key.size = unit(0.04, "npc"),
          legend.spacing.y = unit(.01, "npc"),
          legend.text = element_text(size = rel(1.5)),
          legend.title = element_blank(),
          axis.title = element_text(hjust = 0.5,size = rel(2)),
          axis.text  = element_text(size = rel(1.5))
    )
  # gp_d
  filenm <- paste0("Tunability_", dt_name,"_boxplot.pdf")
  if (save) {
    ggsave(filenm, plot =gp_d, width = width, height = height, device = device)
  }
  
  
  dat_opt <- data.frame(opt = opt_RF_med %>% as.numeric(),
                        mad = opt_RF_mad %>% as.numeric(),
                        dataset = dt_name,
                        Type = c("maxnodes","mtry") %>% as.factor()
  )
  
  
  
  result <- enlist(dt_name, mtry0, maxnodes0, 
                   dat_d, dat_opt, 
                   gp_d
  )
  
  return(result)
  
}






