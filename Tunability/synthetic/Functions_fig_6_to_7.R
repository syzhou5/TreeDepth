
# This file contains functions used in simulations for 
# Figure 6, 7 and B7 in the paper.

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


# The following function is used to aggregate results across
# repeated simulations for each setting and to calculate tunability.

# FinalResult : outputs of simulations from running the wrapper
# g_title     : title of the plot
# mtry0       : default for mtry, as a proprotion of p
# maxnode0    : default for maxnodes, as a proportion of n
# lgd_pos     : position of legends
# ymax_box    : upper limits for y-axis in boxplots
# ymax_scatter: upper limits for y-axis in scatterplots
# filename    : suffix for file names of plots 

summary_and_plot  <- function(FinalResult,
                              g_title,
                              mtry0=0.3, maxnode0=1,
                              lgd_pos = c(0.02, 0.98),
                              ymax_box = NULL,
                              ymax_scatter = NULL,
                              filename = "Tunability",
                              save = T,
                              width = 8,
                              height = 6,
                              device = "pdf"){
  

  
  # Extract simulation parameters 
  temp      <- FinalResult[[1]]
  n         <- temp[['n']]
  nval      <- temp[['nval']]
  p         <- temp[['p']]
  s         <- temp[['s']]
  sigma     <- temp[['sigma']]
  SNR       <- temp[['SNR']]
  l_SNR     <- length(SNR)
  maxnode   <- temp[['maxnode']]
  mtry      <- temp[['mtry']]
  l_maxnd   <- length(maxnode)
  l_mtry    <- length(mtry)
  nrep      <- length(FinalResult)
  ntree     <- temp[['ntree']]
  temp_len  <- length(temp)
  
  temp      <- FinalResult
  # Only keep MSEs
  temp[[1]][1:(temp_len -2)] <- NULL  
  
  
  # Index of the candidate of mtry and maxnode closest to the default value
  mtry_dft  <- which.min(abs(mtry - mtry0))
  depth_dft <- which.min(abs(maxnode - maxnode0*n))
  
  
  # Initialize lists for tunability (d), RTE and optimal values of parameters
  
  d_RF    <- d_mtry <- d_depth <- d_additional <- 
    matrix(NA, nrow = nrep, ncol = l_SNR)
  
  opt_RF  <- matrix(NA, nrow = nrep, ncol = 2*l_SNR)
  colnames(opt_RF) <- rep(c("maxnodes", "mtry"), l_SNR)
  
  colnames(d_RF) <- colnames(d_mtry) <-
    colnames(d_depth) <- colnames(d_additional) <- 
    paste("SNR", SNR)
  
  
  
  for (k in 1:l_SNR) {
    
    for (i in 1:nrep) {
      ## For each repetition at each SNR, 
      ## find the tunability and optimal value of parameters
      
      RTE  <- FinalResult[[i]]$MSE_test[[k]]/sigma[k]^2
      # minimum relative test error, used as reference
      RTE0 <- min(RTE)
      
      # optimal combination of mtry and maxnodes
      ind_opt_RF <- which(RTE == RTE0, arr.ind = T)
      opt_RF[i, (2*k-1):(2*k)] <- c(maxnode[ind_opt_RF[1]]/n, mtry[ind_opt_RF[2]])
      
      d_RF[i,k]    <- RTE[depth_dft, mtry_dft] - RTE0
      d_mtry[i,k]  <- RTE[depth_dft, mtry_dft] - min(RTE[depth_dft, ])
      d_depth[i,k] <- RTE[depth_dft, mtry_dft] - min(RTE[, mtry_dft])
      d_additional[i,k] <- min(c(min(RTE[depth_dft, ]),
                                 min(RTE[, mtry_dft]))) - RTE0
      
    } # end for nrep
    
  } # end for SNR
  
  
  # d_RF_ave is the average of d_RF across 100 repetitions
  d_RF_ave    <- colMeans(d_RF)
  d_RF_se     <- se(d_RF)
  
  d_mtry_ave  <- colMeans(d_mtry)
  d_mtry_se   <- se(d_mtry)
  
  d_depth_ave <- colMeans(d_depth)
  d_depth_se  <- se(d_depth)
  
  d_additional_ave <- colMeans(d_additional)
  d_additional_se  <- se(d_additional)
  
  opt_RF_med  <- apply(opt_RF, 2, median, na.rm=TRUE)%>% matrix(nrow = 2)
  opt_RF_mad   <- (apply(opt_RF, 2, mad, na.rm=TRUE) /
                     sqrt(colSums(!is.na(opt_RF)))) %>% matrix(nrow = 2)
  
  
  rownames(opt_RF_med) <- rownames(opt_RF_mad) <- c("maxnodes","mtry")
  colnames(opt_RF_med) <- colnames(opt_RF_mad) <- paste("SNR", SNR)
  
  
  # Generate boxplots of tunability d
  dat_d <- data.frame(d         = c(d_RF    %>% as.numeric(),
                                    d_mtry  %>% as.numeric(),
                                    d_depth %>% as.numeric(),
                                    d_additional %>% as.numeric()),
                      SNR       = rep(SNR, each = nrep)       %>% as.factor(),
                      Type      = c("RF", "mtry", "maxnodes", "Additional")    %>% rep(each = nrep*l_SNR) %>% as.factor(),
                      mtry_dft  = mtry[mtry_dft] %>% as.factor(),
                      depth_dft = maxnode[depth_dft] %>% as.factor()
  )
  
  
  gp_d <- ggplot(dat_d, aes(x = SNR, y = d)) +
    geom_boxplot(aes(color = Type)) +
    xlab("Signal-to-Noise Ratio") +
    ylab("Tunability") + 
    ggtitle(g_title)
  
  if (!is.null(ymax_box)) {
    gp_d <- gp_d + coord_cartesian(ylim = c(0,ymax_box)) 
  }
  
  gp_d <- gp_d + theme_bw() +
    theme(legend.position = lgd_pos,
          legend.justification = c(0,1),
          plot.title = element_text(hjust = 0.5,size = rel(2.5), face = 'bold'),
          legend.key.size = unit(0.04, "npc"),
          legend.spacing.y = unit(.01, "npc"),
          legend.text = element_text(size = rel(1.5)),
          legend.title = element_blank(),
          axis.title = element_text(hjust = 0.5,size = rel(2)),
          axis.text  = element_text(size = rel(1.5))
    )
  
  filenm <- paste0(filename,"_mtry_", 10*mtry0, "_boxplot.pdf")
  if (save) {
    ggsave(filenm, plot =gp_d, width = width, height = height, device = device)
  }
  
  # Scatterplot of average tunability across 100 repetitions
  dat_d_ave  <- data.frame(d = c(d_RF_ave    %>% as.numeric(),
                                 d_mtry_ave  %>% as.numeric(),
                                 d_depth_ave %>% as.numeric(),
                                 d_additional_ave %>% as.numeric()),
                           se = c(d_RF_se    %>% as.numeric(),
                                  d_mtry_se  %>% as.numeric(),
                                  d_depth_se %>% as.numeric(),
                                  d_additional_se %>% as.numeric()),
                           SNR = SNR ,
                           Type = c("RF", "mtry", "maxnodes", "Additional") %>% rep(each = l_SNR),
                           mtry_dft  = mtry[mtry_dft],
                           depth_dft = maxnode[depth_dft]
  )
  
  gp_d_ave <- ggplot(dat_d_ave, aes(x = SNR, y = d, color = Type)) +
    geom_point() +
    geom_line(lwd = 1) +
    geom_errorbar(aes(ymin = d - se, ymax = d + se)) +
    scale_x_continuous(breaks = SNR,trans = "log") +
    xlab("Signal-to-Noise Ratio") +
    ylab("Tunability") + 
    ggtitle(g_title) 
  
  if (!is.null(ymax_scatter)) {
    gp_d_ave1 <- gp_d_ave1 + coord_cartesian(ylim = c(0, ymax_scatter)) 
  }
  
  
  gp_d_ave <- gp_d_ave + theme_bw() +
    theme(legend.position = lgd_pos,
          legend.justification = c(0,1),
          plot.title = element_text(hjust = 0.5,size = rel(2.5), face = 'bold'),
          legend.key.size = unit(0.04, "npc"),
          legend.spacing.y = unit(.01, "npc"),
          legend.text = element_text(size = rel(1.5)),
          legend.title = element_blank(),
          axis.title = element_text(hjust = 0.5,size = rel(2)),
          axis.text  = element_text(size = rel(1.5))
    )
  # gp_d_ave
  filenm <- paste0(filename,"_mtry_", 10*mtry0, "_scatterplot.pdf")
  
  if (save) {
    ggsave(filenm, plot =gp_d_ave, width = width, height = height, device = "pdf")
  }
  
  dat_opt <- data.frame(opt = opt_RF_med %>% as.numeric(),
                        mad = opt_RF_mad %>% as.numeric(),
                        SNR = rep(SNR, each = 2),
                        Type = c("maxnodes","mtry") %>% as.factor(),
                        mtry_dft  = mtry[mtry_dft],
                        depth_dft = maxnode[depth_dft]
  )
  
  gp_opt <- ggplot(dat_opt, aes(x = SNR, y = opt, color = Type)) + 
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = opt - mad, ymax = opt + mad)) +
    scale_color_manual(values = c("tomato2", "deepskyblue2")) +
    scale_x_continuous(breaks = SNR,trans = "log") +
    xlab("Signal-to-Noise Ratio") +
    scale_y_continuous(
      name = "mtry",
      sec.axis = dup_axis(name = "maxnodes")
    ) +
    coord_cartesian(ylim = c(0,1)) +
    ggtitle(g_title) +
    theme_bw() +
    theme(legend.position = c(0.02, 0.98),
          legend.justification = c(0,1),
          plot.title = element_text(hjust = 0.5,size = rel(2.5), face = 'bold'),
          legend.key.size = unit(0.04, "npc"),
          legend.spacing.y = unit(.01, "npc"),
          legend.text = element_text(size = rel(1.5)),
          legend.title = element_blank(),
          axis.title = element_text(hjust = 0.5,size = rel(2)),
          axis.text  = element_text(size = rel(1.5)),
          axis.title.y.left = element_text(color = "deepskyblue2"),
          axis.title.y.right = element_text(color = "tomato2")
    )
  # gp_opt
  filenm <- paste0(filename,"_OptPara", ".pdf")
  
  
  if (save) {
    ggsave(filenm, plot =gp_opt, width = width, height = height, device = "pdf")
  }
  
  
  result <- enlist(FinalResult,  
                   d_RF, d_RF_ave, d_RF_se,
                   d_mtry, d_mtry_ave, d_mtry_se,
                   d_depth, d_depth_ave, d_depth_se,
                   d_additional, d_additional_ave, d_additional_se,
                   opt_RF, opt_RF_med, opt_RF_mad,
                   mtry0, maxnode0,
                   dat_d, dat_d_ave, dat_opt,
                   gp_d, gp_d_ave, gp_opt
  )
  
  return(result)
}









