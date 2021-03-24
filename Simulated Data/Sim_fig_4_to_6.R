
# This file contains example codes of simulations for 
# Figure 4 to 6.

# Here we consider the case of the low setting with 
# n = 100, p = 10 and s = 5. Others can be obtained 
# by modifying these parameters.


library(randomForest)
library(ggplot2)
library(tidyverse)

# The following is the wrapper function that can be run in
# parallel for repeating the simulation 100 times.

wrapper <- function(nsim){
  
  seed=nsim
  set.seed(seed)
  cat(nsim, "\n")
  
  # Specify parameters for generating data
  # low setting
  n <- 100         
  p <- 10      
  s <- 5         
  rho <- 0.35     
  beta.type <- 2  
  snr <- exp(seq(log(0.05),log(6),length=10))
  nval <- 1000    
  ntest <- 1000
  
  #  Specify parameters for the forest
  nnode <- 10
  ndsize <- c(round(seq(from=n/2,to=5,length.out = nnode)),3,1)
  
  sim.obj  <- sim.rf.optTree(n,p,nval,ntest,rho,s,beta.type,snr,ndsize)
  
  
  if (nsim==1) {
    result <- c(sim.obj,list(n=n,p=p,s=s,rho=rho,beta.type=beta.type,
                             snr=snr,nval=nval,ntest=ntest,ndsize=ndsize))
  }else{
    result <- sim.obj
  }
  
  return(result)
}


# Use the wrapper function above to run the simulation 100 times
obj <- lapply(1:100, wrapper)


n    <- obj[[1]]$n
p    <- obj[[1]]$p
s    <- obj[[1]]$s
nsnr <- length(obj[[1]]$snr)
nrep <- length(obj)

# a matrix of nrep x nsnr, each col corresponding to a SNR
ndsize.tune <- t(sapply(X=obj,FUN=function(x)x[['ndsize.tune']]))
snr         <- round(obj[[1]]$snr,2)


dat <-data.frame(ndsize.tune = as.numeric(ndsize.tune), 
                 SNR         = as.factor(rep(snr,each=nrep)), 
                 n           = rep(n, nrep * nsnr) %>% as.factor(),
                 p           = rep(p, nrep * nsnr) %>% as.factor(),
                 s           = rep(s, nrep * nsnr) %>% as.factor())


gp.box <- ggplot(dat, aes(x=SNR, y=ndsize.tune)) +
  geom_boxplot() +
  scale_x_discrete(breaks = snr, labels = c(0.05, "", 0.14, "", 0.42, "", 1.22, "", 3.52, ""))+
  ggtitle(paste0("n=",n," p=",p," s=",s)) +
  xlab("Signal-to-Noise Ratio") +
  ylab("Optimal Nodesize") +
  geom_smooth(aes(group = 1), se = F) + 
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = rel(2), face = "bold"),
        axis.title = element_text(size = rel(1.7)),
        axis.text = element_text(size = rel(1.5))
  )
gp.box

