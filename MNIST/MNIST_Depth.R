
# This file contains the simulations on the MNIST dataset.
# Functions used below are in the R script Functions_MNIST.R


#### Load Libraries
library(randomForest)
library(magrittr)
library(ggplot2)
library(tidyverse)


load_mnist()
n_tr   <- train$n
n_test <- test$n


# Sample 2000 obs from both training and test sets
n0 <- 2000
set.seed(1)
id_tr    <- sample(1:n_tr, n0, replace = F)
id_test0 <- sample(1:n_test, 2*n0, replace = F)
id_val   <- id_test0[1:n0]
id_test  <- id_test0[(n0+1):(2*n0)]

x_train  <- train$x[id_tr, ]
y_train  <- train$y[id_tr ]

# This validation set is used to choose the maxinodes for building tuned forests
x_val    <- test$x[id_val, ]
y_val    <- test$y[id_val ]
x_test   <- test$x[id_test, ]
y_test   <- test$y[id_test ]

# The job is to predict whether the handwritten digit is a 1
y_out    <- 1

# Respnse for regression
y_train_reg <- ifelse(y_train == y_out, 1, 0)
y_val_reg   <- ifelse(y_val == y_out, 1, 0)
y_test_reg  <- ifelse(y_test == y_out, 1, 0)


# Responses for classification
y_train_class <- factor(y_train_reg)
y_val_class   <- factor(y_val_reg) 
y_test_class  <- factor(y_test_reg)


maxnodes <- c(2:9,
              seq(10, 95, by = 5),
              seq(100, 975, by = 25),
              seq(1000, n0, by = 100))

ntree        <- c(1:4, seq(5, 50, by = 5), 
                  seq(75, 500, by = 25))

# Obtain the training, validation and test error of a randomized tree
# for different jobs. 
# classification/regression + bootstrap samples/original samples

set.seed(1)
tree_class_boot   <- lapply(maxnodes, tree_randomized_error, 
                            replace = T,       job   = "classification", 
                            x       = x_train, y     = y_train_class, 
                            xtest   = x_test,  ytest = y_test_class)

tree_class_noboot <- lapply(maxnodes, tree_randomized_error, 
                            replace = F,       job   = "classification", 
                            x       = x_train, y     = y_train_class, 
                            xtest   = x_test,  ytest = y_test_class)

tree_reg_boot     <- lapply(maxnodes, tree_randomized_error, 
                            replace = T,       job   = "regression", 
                            x       = x_train, y     = y_train_reg, 
                            xtest   = x_test,  ytest = y_test_reg)

tree_reg_noboot   <- lapply(maxnodes, tree_randomized_error, 
                            replace = F,       job   = "regression", 
                            x       = x_train, y     = y_train_reg, 
                            xtest   = x_test,  ytest = y_test_reg)

tree_out <- list(tree_class_boot, tree_class_noboot, 
                 tree_reg_boot,   tree_reg_noboot)


tree_error <- lapply(tree_out, 
                     function(x){
                       error_train_tree <- sapply(x, function(y){y[1]})
                       error_test_tree  <- sapply(x, function(y){y[2]})
                       enlist(error_train_tree, error_test_tree)
                     })


tree_class_boot_val   <- lapply(maxnodes, tree_randomized_error, 
                                replace = T,       job   = "classification", 
                                x       = x_train, y     = y_train_class, 
                                xtest   = x_val,   ytest = y_val_class)

tree_class_noboot_val <- lapply(maxnodes, tree_randomized_error, 
                                replace = F,       job   = "classification", 
                                x       = x_train, y     = y_train_class, 
                                xtest   = x_val,   ytest = y_val_class)

tree_reg_boot_val     <- lapply(maxnodes, tree_randomized_error, 
                                replace = T,       job   = "regression", 
                                x       = x_train, y     = y_train_reg, 
                                xtest   = x_val,   ytest = y_val_reg)

tree_reg_noboot_val   <- lapply(maxnodes, tree_randomized_error, 
                                replace = F,       job   = "regression", 
                                x       = x_train, y     = y_train_reg, 
                                xtest   = x_val,   ytest = y_val_reg)

tree_out_val <- list(tree_class_boot_val, tree_class_noboot_val, 
                     tree_reg_boot_val,   tree_reg_noboot_val)


tree_error_val <- lapply(tree_out_val, 
                         function(x){
                           error_train_tree <- sapply(x, function(y){y[1]})
                           error_test_tree  <- sapply(x, function(y){y[2]})
                           enlist(error_train_tree, error_test_tree)
                         })



# Optimal maxnodes for a randomized tree based on the validation set above.

maxnodes_opt <- lapply(tree_error_val,
                       function(x){
                         max_opt <- maxnodes[which.min(x[[2]]) ]
                       })
names(tree_error_val) <- 
  names(tree_error)  <- names(maxnodes_opt) <- c("Class + Boot",
                                                 "Class + NoBoot",
                                                 "Reg + Boot",
                                                 "Reg + NoBoot")





# Obtain the training and test error for full-depth forests, tuned forests, 
# and shallow forests.
# The following is for classification forests with bootstrap samples.
# Others can be done similarly.

forests <- lapply(1:length(ntree), 
                  function(nsim){
                    
                    set.seed(maxnodes_opt$`Class + Boot` - 1 + nsim)
                    
                    # tuned forests
                    rf      <- randomForest(x = x_train,
                                            y = y_train_class,
                                            xtest = x_test,
                                            ytest = y_test_class,
                                            ntree = ntree[nsim],
                                            maxnodes = maxnodes_opt$`Class + Boot`,
                                            nodesize = 1,
                                            keep.forest = T) 
                    
                    pred_train  <- predict(rf, newdata = x_train)
                    error_train <- mean(pred_train != y_train_class  )
                    
                    pred_test   <- predict(rf, newdata = x_test)
                    error_test  <- mean(pred_test != y_test_class )
                    
                    # full-depth forests
                    rf_full <- randomForest(x = x_train,
                                            y = y_train_class,
                                            xtest = x_test,
                                            ytest = y_test_class,
                                            ntree = ntree[nsim],
                                            nodesize = 1,
                                            keep.forest = T)
                    
                    pred_full_train <- predict(rf_full, newdata = x_train)
                    error_full_train <- mean(pred_full_train != y_train_class  )
                    
                    pred_full_test <- predict(rf_full, newdata = x_test)
                    error_full_test<- mean(pred_full_test != y_test_class )
                    
                    # shallow random forest
                    rf_s    <- randomForest(x = x_train,
                                            y = y_train_class,
                                            xtest = x_test,
                                            ytest = y_test_class,
                                            ntree = ntree[nsim],
                                            maxnodes = 10,
                                            nodesize = 1,
                                            keep.forest = T) 
                    
                    pred_s_train  <- predict(rf_s, newdata = x_train)
                    error_s_train <- mean(pred_s_train != y_train_class  )
                    
                    pred_s_test   <- predict(rf_s, newdata = x_test)
                    error_s_test  <- mean(pred_s_test != y_test_class )
                    
                    
                    result <- enlist(error_full_test, error_full_train, 
                                     error_test,      error_train,
                                     error_s_test,    error_s_train)
                    
                    
                    return(result)
                  })





# Smooth error curves of a randomized tree
tree_error_lowess <- lapply(tree_error,
                            function(x){
                              temp <- lapply(x, function(a){
                                a_lowess <- lowess(x = maxnodes, y =a, f = 0.3)
                                a_out <- a_lowess$y
                                return(a_out)
                              })
                            })


date <- " "
job <- "class"
Boot <- T
lgd <- F


rf_full_error <- rf_tuned_error <- rf_shallow_error <- 
  dat_tree <- dat_rf <- dat <- gp <- list()

plot_depth(forests, date, job, Boot, lgd)







