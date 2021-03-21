
# This file contains functions used in the simulations on the MNIST dataset.


# load_mnist and show_digit are from https://gist.github.com/brendano/39760
# for loading the MNIST dataset into R
# The MNIST dataset is downloaded from http://yann.lecun.com/exdb/mnist/
load_mnist <- function() {
  load_image_file <- function(filename) {
    ret = list()
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    ret$n = readBin(f,'integer',n=1,size=4,endian='big')
    nrow = readBin(f,'integer',n=1,size=4,endian='big')
    ncol = readBin(f,'integer',n=1,size=4,endian='big')
    x = readBin(f,'integer',n=ret$n*nrow*ncol,size=1,signed=F)
    ret$x = matrix(x, ncol=nrow*ncol, byrow=T)
    close(f)
    ret
  }
  load_label_file <- function(filename) {
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    n = readBin(f,'integer',n=1,size=4,endian='big')
    y = readBin(f,'integer',n=n,size=1,signed=F)
    close(f)
    y
  }
  train <<- load_image_file('train-images.idx3-ubyte')
  test <<- load_image_file('t10k-images.idx3-ubyte')
  
  train$y <<- load_label_file('train-labels.idx1-ubyte')
  test$y <<- load_label_file('t10k-labels.idx1-ubyte')  
}


show_digit <- function(arr784, col=gray(12:1/12), ...) {
  image(matrix(arr784, nrow=28)[,28:1], col=col, ...)
}



# The following function enlist is from the bestsubset package.
# https://github.com/ryantibs/best-subset/
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



# The following function is for calculating the training 
# and test error of a randomized tree.
tree_randomized_error <- function(maxnodes, x, y, xtest, 
                                  ytest, replace = T, 
                                  job = "classification"){
  
  if (!job %in% c("classification", "regression")) {
    stop("Job must be either classification or regression")
  }
  
  if (job == "classification") {
    if (replace) {
      # Using bootstrap samples
      tree <- randomForest(x        = x,
                           y        = y,
                           xtest    = xtest,
                           ytest    = ytest,
                           maxnodes = maxnodes,
                           ntree    = 1,
                           nodesize = 1,
                           keep.forest = T)
    }else{
      # Using original samples
      tree <- randomForest(x        = x,
                           y        = y,
                           xtest    = xtest,
                           ytest    = ytest,
                           maxnodes = maxnodes,
                           replace  = F,
                           sampsize = n0,
                           ntree    = 1,
                           nodesize = 1, 
                           keep.forest = T)
    }
    pred_train <- predict(tree, newdata = x)
    error_train <- mean(pred_train != y)
    
    pred_test <- predict(tree, newdata = xtest)
    error_test<- mean(pred_test != ytest)
    
    
  }else{
    
    if (replace) {
      # Using bootstrap samples
      tree <- randomForest(x        = x,
                           y        = y,
                           xtest    = xtest,
                           ytest    = ytest,
                           maxnodes = maxnodes,
                           ntree    = 1,
                           nodesize = 1, 
                           keep.forest = T)
    }else{
      # Using original samples
      tree <- randomForest(x        = x,
                           y        = y,
                           xtest    = xtest,
                           ytest    = ytest,
                           maxnodes = maxnodes,
                           replace  = F,
                           sampsize = n0,
                           ntree    = 1,
                           nodesize = 1, 
                           keep.forest = T)
    }
    
    pred_train <- predict(tree, newdata = x)
    error_train <- mean((pred_train - y)^2 )
    
    pred_test <- predict(tree, newdata = xtest)
    error_test<- mean((pred_test - ytest)^2 )
    
  }
  
  
  error <- c(error_train, error_test)
  return(error)
}



# The function for generating plots of error against model complexity

plot_depth <- function(forests, date, job = c("class", "reg"), Boot = T, lgd = T){
  
  if (job == "class" & Boot == T) {
    name <- "Class + Boot"
    y_lab <- "0-1 Loss"
    plot_title <- "Classification with Bootstrap Samples"
    filename <- paste0(date, "_Class_Boot")
    
  } else if (job =="class" & Boot ==F){
    name <- "Class + NoBoot"
    y_lab <- "0-1 Loss"
    plot_title <- "Classification with Original Samples"
    filename <- paste0(date, "_Class_NoBoot")
  } else if (job == "reg" & Boot == T){
    name <- "Reg + Boot"
    y_lab <- "L2 Loss"
    plot_title <- "Regression with Bootstrap Samples"
    filename <- paste0(date, "_Reg_Boot")
  }else{
    name <- "Reg + NoBoot"
    y_lab <- "L2 Loss"
    plot_title <- "Regression with Original Samples"
    filename <- paste0(date, "_Reg_NoBoot")
  }
  
  
  rf_full_error[[name]] <- list()
  
  rf_full_error[[name]][["error_train"]]  <- sapply(forests, 
                                                    function(x){
                                                      x$error_full_train
                                                    })
  rf_full_error[[name]][["error_test"]]   <- sapply(forests, 
                                                    function(x){
                                                      x$error_full_test
                                                    })
  
  rf_tuned_error[[name]] <- list()
  rf_tuned_error[[name]][["error_train"]] <- sapply(forests, 
                                                    function(x){
                                                      x$error_train
                                                    })
  rf_tuned_error[[name]][["error_test"]]  <- sapply(forests, 
                                                    function(x){
                                                      x$error_test
                                                    })
  
  rf_shallow_error[[name]] <- list()
  rf_shallow_error[[name]][["error_train"]] <- sapply(forests, 
                                                      function(x){
                                                        x$error_s_train
                                                      })
  rf_shallow_error[[name]][["error_test"]]  <- sapply(forests, 
                                                      function(x){
                                                        x$error_s_test
                                                      })
  
  dat_tree[[name]]    <- data.frame(maxnodes = maxnodes %>% rep(2*2),
                                    Error = c(tree_error[[name]] %>% unlist(),
                                              tree_error_lowess[[name]] %>% unlist()),
                                    Type = (c("Train", "Test") %>% rep(each = length(maxnodes))) %>% rep(2) %>% factor(),
                                    ToF = (c("Tree") %>% rep(each = 2*length(maxnodes))) %>% rep(2) %>% factor(),
                                    group = rep(1:4, each =  length(maxnodes)),
                                    Alpha = rep(c(0.95, 1), each = 2* length(maxnodes)),
                                    LOWESS = c("N", "Y") %>% rep(each = 2* length(maxnodes)),
                                    x_lab = (1: length(maxnodes)) %>% rep(4)) 

  dat_rf[[name]]    <- data.frame(Error = c(rf_tuned_error[[name]] %>% unlist(),
                                            rf_full_error[[name]] %>% unlist(),
                                            rf_shallow_error[[name]] %>% unlist()),
                                  Type = c(rep("Train", length(ntree)),
                                           rep("Test", length(ntree))) %>% rep(3) %>% factor(),
                                  ToF = c(rep("RF Tuned", 2*length(ntree)),
                                          rep("RF Full", 2*length(ntree)),
                                          rep("RF Shallow", 2*length(ntree))) %>% factor(),
                                  ntree = rep(ntree, 6),
                                  group = rep(5:10, each = length(ntree)),
                                  Alpha = rep(1, 6*length(ntree)),
                                  LOWESS = rep("N", 6*length(ntree)),
                                  x_lab = c((rep(1:length(ntree), 2))*2 + which(maxnodes == maxnodes_opt[[name]]) - 2,
                                            (rep(1:length(ntree), 2))*2 + length(maxnodes) - 2,
                                            (rep(1:length(ntree), 2))*2 + which(maxnodes == 10) - 2)
  )
  
  dat[[name]]       <- bind_rows(dat_tree[[name]] , dat_rf[[name]] )
  
  # Full-depth Forests
  dat_full <- subset(dat[[name]], ToF %in% c("Tree", "RF Full"))
  gp_full  <- ggplot(dat_full, aes(x = x_lab, y = Error, color = ToF)) + 
    geom_line(aes(group = group, linetype = Type, alpha = Alpha)) +
    geom_vline(xintercept = c(dat_full$x_lab[which(dat_full$maxnodes == 10)][1],
                              dat_full$x_lab[which(dat_full$maxnodes == 100)][1],
                              dat_full$x_lab[which(dat_full$maxnodes == 1000)][1],
                              dat_full$x_lab[which(dat_full$maxnodes == 2000)][1],
                              dat_full$x_lab[which(dat_full$maxnodes == maxnodes_opt[[name]])][1],
                              dat_full$x_lab[which(dat_full$maxnodes == 10)][1],
                              dat_full$x_lab[which(dat_full$maxnodes == 2000)][1]),
               linetype   = c(rep(3, 6), 5), 
               color      = c(rep("black", 4), "deepskyblue", "olivedrab", "red"),
               size       = c(rep(0.5, 6), rep(1, 1))) +
    scale_alpha_continuous(range = c(0.3,1), guide = F) +
    scale_color_manual(values = c("purple", "red")) +
    scale_x_continuous(breaks = c(dat_full$x_lab[which(dat_full$maxnodes == 2)][1],
                                  dat_full$x_lab[which(dat_full$maxnodes == 50)][1],
                                  dat_full$x_lab[which(dat_full$maxnodes == 500)][1],
                                  dat_full$x_lab[which(dat_full$maxnodes == 2000)][1],
                                  dat_full$x_lab[which(dat_full$ntree == 50)][1],
                                  dat_full$x_lab[which(dat_full$ntree == 400)][1]),
                       labels = c("2/1", "50/1", "500/1", "2000/1", 
                                  "2000/50", "2000/400")) +
    xlab("Maxnodes/Trees") +
    ylab(y_lab) + 
    ggtitle(plot_title) +
    ylim(dat[[name]]$Error %>% min() - 0.01,
         dat[[name]]$Error %>% max() + 0.01) +
    theme_bw() +
    theme(legend.position = c(0.98, 0.98),
          legend.justification = c(1,1),
          plot.title = element_text(hjust = 0.5,size = rel(2.3), face = 'bold'),
          legend.key.size = unit(0.04, "npc"),
          legend.spacing.y = unit(.01, "npc"),
          legend.text = element_text(size = rel(2)),
          legend.title = element_text(size = rel(2)),
          axis.title = element_text(hjust = 0.5,size = rel(2)),
          axis.text  = element_text(size = rel(1.5))
    )
  if (lgd == F) {
    gp_full <- gp_full + theme(legend.position = "none")
  }
  
  ggsave(filename = paste0(filename, "_RFfull.pdf"), plot = gp_full, device = "pdf", width = 6, height = 6)
  
  # Tuned Forests
  dat_tuned <- subset(dat[[name]], (ToF =="Tree" & maxnodes <=  maxnodes_opt[[name]]) | ToF == "RF Tuned")
  gp_tuned  <- ggplot(dat_tuned, aes(x = x_lab, y = Error, color = ToF)) + 
    geom_line(aes(group = group, linetype = Type, alpha = Alpha)) +
    geom_vline(xintercept = dat_tuned$x_lab[which(dat_tuned$maxnodes == maxnodes_opt[[name]])][1],
               linetype   = 5, 
               color      = "deepskyblue",
               size       = 1) +
    scale_alpha_continuous(range = c(0.3,1), guide = F) +
    scale_color_manual(values = c("purple", "deepskyblue")) +
    scale_x_continuous(breaks = c(dat_tuned$x_lab[which(dat_tuned$maxnodes == 2)][1],
                                  dat_tuned$x_lab[which(dat_tuned$maxnodes == maxnodes_opt[[name]])][1],
                                  dat_tuned$x_lab[which(dat_tuned$ntree == 50)][1],
                                  dat_tuned$x_lab[which(dat_tuned$ntree == 400)][1]),
                       labels = c("2/1", 
                                  paste0(maxnodes_opt[[name]], "/1"),
                                  paste0(maxnodes_opt[[name]], "/50"),
                                  paste0(maxnodes_opt[[name]], "/400"))) +
    xlab("Maxnodes/Trees") +
    ylab(y_lab) + 
    ggtitle(plot_title) +
    ylim(dat[[name]]$Error %>% min() - 0.01,
         dat[[name]]$Error %>% max() + 0.01) +
    theme_bw() +
    theme(legend.position = c(0.98, 0.98),
          legend.justification = c(1,1),
          plot.title = element_text(hjust = 0.5,size = rel(2.3), face = 'bold'),
          legend.key.size = unit(0.04, "npc"),
          legend.spacing.y = unit(.01, "npc"),
          legend.text = element_text(size = rel(2)),
          legend.title = element_text(size = rel(2)),
          axis.title.x = element_text(hjust = 0.5,size = rel(2)),
          axis.title.y = element_blank(),
          axis.text  = element_text(size = rel(1.5))
    )
  if (lgd == F) {
    gp_tuned <- gp_tuned + theme(legend.position = "none")
  }
  ggsave(filename = paste0(filename, "_RFtuned.pdf"), plot = gp_tuned, device = "pdf", width = 6, height = 6)
  
  # Shallow Forests
  dat_shallow <- subset(dat[[name]], (ToF =="Tree" & maxnodes <=  10) | ToF == "RF Shallow")
  gp_shallow  <- ggplot(dat_shallow, aes(x = x_lab, y = Error, color = ToF)) + 
    geom_line(aes(group = group, linetype = Type, alpha = Alpha)) +
    geom_vline(xintercept = dat_shallow$x_lab[which(dat_shallow$maxnodes == 10)][1],
               linetype   = 5, 
               color      = "olivedrab",
               size       = 1) +
    scale_alpha_continuous(range = c(0.3,1), guide = F) +
    scale_color_manual(values = c("purple", "olivedrab")) +
    scale_x_continuous(breaks = c(dat_shallow$x_lab[which(dat_shallow$maxnodes == 2)][1],
                                  dat_shallow$x_lab[which(dat_shallow$maxnodes == 10)][1],
                                  dat_shallow$x_lab[which(dat_shallow$ntree == 50)][1],
                                  dat_shallow$x_lab[which(dat_shallow$ntree == 400)][1]),
                       labels = c("2/1", 
                                  "10/1", 
                                  "10/50",
                                  "10/400")) +
    xlab("Maxnodes/Trees") +
    ylab(y_lab) + 
    ggtitle(plot_title) +
    ylim(dat[[name]]$Error %>% min() - 0.01,
         dat[[name]]$Error %>% max() + 0.01) +
    theme_bw() +
    theme(legend.position = c(0.98, 0.98),
          legend.justification = c(1,1),
          plot.title = element_text(hjust = 0.5,size = rel(2.3), face = 'bold'),
          legend.key.size = unit(0.04, "npc"),
          legend.spacing.y = unit(.01, "npc"),
          legend.text = element_text(size = rel(2)),
          legend.title = element_text(size = rel(2)),
          axis.title.x = element_text(hjust = 0.5,size = rel(2)),
          axis.title.y = element_blank(),
          axis.text  = element_text(size = rel(1.5))
    )
  
  if (lgd == F) {
    gp_shallow <- gp_shallow + theme(legend.position = "none")
  }
  ggsave(filename = paste0(filename, "_RFshallow.pdf"), plot = gp_shallow, device = "pdf", width = 6, height = 6)
  
  
}














