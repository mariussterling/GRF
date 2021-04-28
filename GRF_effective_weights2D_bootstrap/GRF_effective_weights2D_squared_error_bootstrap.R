#install.packages("matrixStats")
rm(list=ls())
library(grf)
library(ggplot2)
library(glue)
library(rstudioapi)
library(parallel)
library(matrixStats)
set.seed(42)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

# Function definition for data --------------------------------------------
get_x = function(n,seed){
  set.seed(seed)
  X= expand.grid(X1 = seq(-0.5, 0.5, length.out = n), X2 = seq(0, 1, 0.1))
  return(X)
}

theta_triangle = function(x, width){
  pmax(1 - abs((x[[1]]) / width), 0)
} 

get_y = function(X, theta, sigma, seed=NULL, reps = 1){
  set.seed(NULL)
  n = nrow(X)
  return(replicate(reps,(theta(X) + rnorm(n, 0, sigma))))
} 

# Data initializing  ------------------------------------------------------
tau = c(0.5)
sig = 1
width = 0.2
c = 1
b = 100 #number of bootstraps for MBS
reps = 20 #repititions for confidence interval

rfs = list()
T_stats = list()
CIs = list()
CIs_std = list()

# Computation ________________________________
for ( n in c(50)){
  T_stat = list()
  set.seed(100)
  ptm <- proc.time()
  X = get_x(n, seed=42) #no replications for X because X is deterministic
  theta_fun = function(X) theta_triangle(X, width)
  theta_true = theta_triangle(X, width)
  Y = get_y(X, theta_fun, sig, NULL,reps)
  rand_for =  function(j)  grf::quantile_forest( X ,data.matrix(Y[,j]), quantile = 0.5)
  rf = lapply(1:reps, rand_for)
  w = lapply(1:reps,function(j) get_sample_weights(rf[[j]]))
  alpha = lapply(1:reps, function(j) w[[j]][1,]) # extracting 1st column of ws as vector for 
  objective_fun = function(theta_opt,Y,alpha) (1/2) * sum( unlist(alpha) * (Y-theta_opt)^(2))
  theta_hat = lapply(1:reps, function(j) sapply(1:nrow(X), function(k)
    optimize(f=objective_fun,interval = c(0,1),
             tol = 0.0001, Y=Y[,j], alpha=w[[j]][k,])[1]))
  
  #kde = lapply(1:reps, function (j) density(Y[,j], n=n)) #estimation of the density
  #f_Y= sapply(1:reps, function(j) unlist(approx(kde[[j]][["x"]], kde[[j]][["y"]], xout = c(theta_hat[[1]]))[2]))
  #psi_fun = function(Y, theta) Y-theta
  #psi = psi_fun(Y,theta_hat)   
  #V_hat = 1/f_Y
  V_hat = -1
  #H_hat = sapply(1:reps, function(j)  n *((var(alpha[[j]] * (Y[,j]-unlist(theta_hat[j]))))))
  H_hat = sapply(1:reps, function(j)  var((w[[j]] %*% (Y[,j]-unlist(theta_hat[j])))@x))
  
  sigma_hat = (H_hat)^(1/2)
  e_multipliers = lapply(1:reps, function(j) lapply(1:b, function(j) rnorm(n, 0, 1)))
  
  #   T_stat = sapply(1:reps, function(k) sapply(1:b, function(j) 
  #    (H_hat[k]^(-1/2)) * sum( (Y[,k]- unlist(theta_hat[k])) * alpha[[k]] * e_multipliers[[k]][[j]])))
  
  T_stat = lapply(1:reps, function(k) sapply(1:b, function(j) 
    (w[[k]] %*% ((H_hat[k]^(-1/2)) * (Y[,k]- unlist(theta_hat[k]))  * e_multipliers[[k]][[j]]))@x))
  
  ## Confidence interval with test stat
  alpha_sig = 0.05
  T_stat_abs = lapply(1:reps, function(j) abs(t(T_stat[[j]])))
  q_star = lapply(1:reps, function(j) colQuantiles(T_stat_abs[[j]] , probs= c(1-alpha_sig)))
  #CI = list(unlist(theta_hat)-q_star*sigma_hat, unlist(theta_hat)+q_star*sigma_hat)
  CI = lapply(1:reps, function(j) list(unlist(theta_hat[[j]])-(q_star[[j]]*sigma_hat[j]),
                                       unlist(theta_hat[[j]])+(q_star[[j]]*sigma_hat[j])))
  
  ## Confidence interval with standard normal
  #q_norm = qnorm(1-alpha_sig) 
  #CI_std = list(unlist(theta_hat)-q_norm*sigma_hat, unlist(theta_hat)+q_norm*sigma_hat)
  #CI_std = do.call("cbind",CI_std)
  #print(proc.time() - ptm)
}

#T_stats[[as.character(n)]]  = T_stat
CIs[[as.character(n)]]  = CI


# Plotting _______________________________________

## plot with multiple confidence intervals 
for (num in c(1,4)){
  for (x2 in c(0.3,0.5)){
    png(file = glue('CI_','n{formatC(as.integer(n), width=4, flag="0")}_',
                    'reps{formatC(as.integer(num) ,flag="0")}_',
                    'x2{formatC(x2*10, width=3 ,flag="0")}_',
                    '.png'),
        width=1500, height=1500)
    
    plot(1, type="n", xlab="X", ylab=bquote(theta), xlim=c(-0.5, 0.5), ylim= c(-0.5,1.25))
    for (i in 1:num){
      pd = data.frame(X1=X$X1,X2=X$X2, theta_hat =  unlist(theta_hat[[i]]),
                      theta_true = theta_true, CI_L = CI[[i]][[1]], CI_U = CI[[i]][[2]] )
      pd = pd[pd$X2==sprintf("%0.1f", x2),]
      pd = pd[order(pd$X1),]
      points(pd$X1, pd$theta_true ,
             ylim=range(pd$theta_true,pd$theta_hat, pd$CI_L, pd$CI_U),
             col='red', main="Confidence intervals", pch=19, cex=2)
      points(pd$X1, pd$theta_hat, col='blue', pch=19, cex=2)
      lines(pd$X1, pd$CI_L,      col='black',pch = 19,  type = "b", lty = 2, cex=0.8)
      lines(pd$X1, pd$CI_U,      col='black',pch = 19,  type = "b", lty = 2, cex=0.8)
    }
    dev.off()
  }}


## plot with average confidence intervals
avg_ci_lower = 0
avg_ci_upper = 0
avg_theta_hat = 0
for (i in seq(1, reps,1)){
  avg_ci_lower = avg_ci_lower + CI[[i]][[1]] 
  avg_ci_upper = avg_ci_upper + CI[[i]][[2]]
  avg_theta_hat = avg_theta_hat + unlist(theta_hat[[i]])
}
avg_ci_lower = avg_ci_lower/reps
avg_ci_upper = avg_ci_upper/reps
avg_theta_hat = avg_theta_hat/reps

for (x2 in c(0.3,0.5)){
  png(file = glue('CI_averaged_','n{formatC(as.integer(n), width=4, flag="0")}_',
                  'reps{formatC(as.integer(reps) ,flag="0")}_',
                  'x2{formatC(x2*10, width=3 ,flag="0")}_',
                  '.png'),
      width=1500, height=1500)
  
  plot(1, type="n", xlab="X", ylab=bquote(theta), xlim=c(-0.5, 0.5), ylim= c(-0.5,1.25))
  pd = data.frame(X1=X$X1,X2=X$X2, theta_hat =  avg_theta_hat,
                  theta_true = theta_true, CI_L = avg_ci_lower, CI_U = avg_ci_upper)
  pd = pd[pd$X2==sprintf("%0.1f", x2),]
  pd = pd[order(pd$X1),]
  points(pd$X1, pd$theta_true ,
         ylim=range(pd$theta_true,pd$theta_hat, pd$CI_L, pd$CI_U),
         col='red', main="Confidence intervals", pch=19, cex=2)
  points(pd$X1, pd$theta_hat, col='blue', pch=19, cex=2)
  lines(pd$X1, pd$CI_L,      col='black',pch = 19,  type = "b", lty = 2, cex=0.8)
  lines(pd$X1, pd$CI_U,      col='black',pch = 19,  type = "b", lty = 2, cex=0.8)
  dev.off()
}
