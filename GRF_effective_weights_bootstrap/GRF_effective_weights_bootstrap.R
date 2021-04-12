rm(list=ls())
library(grf)
library(ggplot2)
library(glue)
library(rstudioapi)
library(parallel)
set.seed(42)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

# Function definition for data --------------------------------------------
get_x = function(n, c, reps =1){
  X = matrix(sort(rep(c, n)), nrow = n, ncol = reps)
  return(X)
}

get_y = function(n, theta, sigma, seed=NULL, reps = 1){
  set.seed(NULL)
  return(replicate(reps,theta + rnorm(n, 0, sigma)))
}
V_func  = function(x, mu, sig, clip = NULL){
  res = -dnorm(x, mean = mu, sd = sig)
  if (!is.null(clip))
    res = sign(res) * pmax(abs(res), clip)
  return(res)
}
'#eps_tilde = function(X, Y, theta, sig, tau, clipV=NULL){
  res = matrix(0, nrow = nrow(X), ncol = length(tau))
  V = - V_func(theta(X), theta(X), sig, clip=clipV) ^ (-1)
  I = Y <= theta(X)
  for (i in 1:length(tau)){
    res[, i] = V * (tau[i] - I)
  }
  return(res)
}'

tau = c(0.5)
sig = 1
theta =  0.9
c = 1
b = 100 #number of bootstraps for MBS
reps = 100 #repititions for confidence interval

rfs = list()
T_stats = list()
CIs = list()
CIs_std = list()
for ( n in c(500,1000)){
  T_stat = list()
  i = 1
  set.seed(100)
  ptm <- proc.time()
  X = get_x(n, c, reps)
  Y = get_y(n, theta, sig, 42, reps)
  rand_for =  function(j)  grf::quantile_forest(data.matrix(X[,j]),data.matrix(Y[,j]), quantile = 0.5)
  rf = lapply(1:ncol(X), rand_for)
  w = lapply(1:ncol(X),function(x) get_sample_weights(rf[[x]]))
  alpha = lapply(1:ncol(X), function(y) w[[y]]@x[1:n]) # extracting 1st column of ws as vector for 
  objective_fun = function(theta,Y,alpha,tau)  sum(((Y-theta)) * as.matrix(alpha) * (tau -   (Y <= theta)))
  theta_hat = sapply(1:ncol(X), function(j) optimize(f=objective_fun,interval = c(0,1), tol = 0.0001, Y=Y[,j], alpha=alpha[[j]], tau=tau)[1])
  
  kde = lapply(1:ncol(X), function (j) density(Y[,j], n=n)) #estimation of the density
  f_Y= sapply(1:ncol(X), function(j) unlist(approx(kde[[j]][["x"]], kde[[j]][["y"]], xout = c(theta_hat[j]))[2]))
  V_hat = 1/f_Y
  H_hat = sapply(1:ncol(X), function(j)  n *((var(alpha[[j]] * (tau -   (Y[,j] <= theta_hat[j]))))))
  sigma_hat = (f_Y^(-2)*H_hat)^(1/2)
  e_multipliers = lapply(1:reps, function(j) lapply(1:b, function(j) rnorm(n, 0, 1)))
  
  T_stat = sapply(1:ncol(X), function(k) sapply(1:b, function(j) -(H_hat[k]^(-1/2)) * sum((tau-(Y[,k] <= theta_hat[k])) * alpha[[k]] * e_multipliers[[k]][[j]])))
  
  ## Confidence interval with test stat
  alpha_sig = 0.05
  T_stat_abs = sapply(1:ncol(X), function(j) abs((T_stat[j,])))
  q_star = sapply(1:reps, function(j) quantile((T_stat_abs[j,]) , 1-alpha_sig))
  CI = list(unlist(theta_hat)-q_star*sigma_hat, unlist(theta_hat)+q_star*sigma_hat)
  CI = do.call("cbind",CI)
  
  ## Confidence interval with standard normal
  q_norm = qnorm(1-alpha_sig) 
  CI_std = list(unlist(theta_hat)-q_norm*sigma_hat, unlist(theta_hat)+q_norm*sigma_hat)
  CI_std = do.call("cbind",CI_std)
  print(proc.time() - ptm)
  # counting the number of times intervals cover the theta
  count = 0 #counting coverage of CI
  count_std = 0 #counting covergae of std normal CI
  for (i in c(seq(1,reps,1))){
    if ((theta > CI[i,1]) & (theta < CI[i,2])  )
      count = count + 1
    if ((theta > CI_std[i,1]) & (theta < CI_std[i,2])  )
      count_std = count_std + 1
  }
  print(count)
  print(count_std)
  
  T_stats[[as.character(n)]]  = T_stat
  CIs[[as.character(n)]]  = CI
  CIs_std[[as.character(n)]]  = CI_std
  d = density(unlist(T_stat[1,]), n=b)
  
  fn = glue(
    'qRF_location_',
    'q{formatC(tau*100, width=3, flag="0")}___',
    'new_var',
    'sigma{formatC(sig*100, width=3, flag="0")}___',
    'theta{formatC(theta*100, width=3, flag="0")}___',
    'n{formatC(as.integer(n), width=4, flag="0")}___',
    'bs{formatC(as.integer(b), width=4, flag="0")}',
    '.png'
  )
  
  png(file = fn)
  
  plot(d, ylab='', xlab='',ylim = c(0,0.45), main="Density plot of bootsrap test statistic")
  ylim = c(0,1)
  x_std = rnorm(n, mean = 0, sd= 1)
  std_d = density(x_std, n=b)
  lines(std_d, col = 'blue')
  dev.off()
  
  png(file = glue('CI_','n{formatC(as.integer(n), width=4, flag="0")}_','.png'))
  
  plot(rep(theta,reps), type = 'l', ylab='', xlab='', ylim=range(theta, CI[,1], CI[,2]), col='black', main="Confidence intervals")
  lines(as.matrix(theta_hat), col='green')
  lines(CI[,1], col='red', lty=2)
  lines(CI[,2], col='red', lty=2)
  dev.off()
  
  png(file = glue('CI_std_','n{formatC(as.integer(n), width=4, flag="0")}_','.png'))
  
  plot(rep(theta,reps), type = 'l', ylab='', xlab='', ylim=range(theta, CI[,1], CI[,2]), col='black', main="Std Confidence intervals")
  lines(as.matrix(theta_hat), col='green')
  lines(CI_std[,1], col='red', lty=2)
  lines(CI_std[,2], col='red', lty=2)
  dev.off()
}

