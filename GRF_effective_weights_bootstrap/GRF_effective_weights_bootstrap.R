rm(list=ls())
library(grf)
library(ggplot2)
library(glue)
library(rstudioapi)
set.seed(42)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

# Function definition for data --------------------------------------------
get_x = function(n, c){
  X = matrix(sort(rep(c, n)), nrow = n)
  return(X)
}

get_y = function(n, theta, sigma, seed=NULL){
  set.seed(NULL)
  return(theta + rnorm(n, 0, sigma))
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
bootstraps = 100

rfs = list()
T_stats = list()
CIs = list()
for ( n in c(500,1000)){
  T_stat = list()
  i = 1
  set.seed(100)
  for (i in c(seq(1,bootstraps,1))){
    X = get_x(n, c)
    Y = get_y(n, theta, sig, 42)
    rfs[[as.character(n)]] = rf = grf::quantile_forest(X, Y, quantiles = tau)
    w = get_sample_weights(rf)
    alpha = as.vector(w[1,]) # extracting 1st column of ws as vector
    ##theta_tilde = theta(X) + alpha %*% eps_tilde(X, Y, theta, sig, tau)
    objective_fun = function(theta,Y,alpha,tau)  sum(((Y-theta)) * as.matrix(alpha) * (tau -   (Y <= theta)))
    theta_hat = optimize(f=objective_fun,interval = c(0,1), tol = 0.0001, Y=Y, alpha=alpha, tau=tau)
    theta_hat = unlist(theta_hat[1])
    kde = density(Y, n=n) #estimation of the density
    f_Y= unlist(approx(kde$x, kde$y, xout = c(theta_hat))[2]) #evaluation of density at theta_hat
    V_hat = - f_Y^(-1) #V funtion
    H_hat = sum(alpha^2 * var((tau -   (Y <= theta_hat))))
    H_hat_1 = n *((var(alpha * (tau -   (Y <= theta_hat))))) #alternative way of calculating H_hat
    e_multipliers = rnorm(n, 0, 1)
    T_star = -(H_hat^(-1/2)) * sum((tau-(Y <= theta_hat)) * alpha * e_multipliers)
    T_stat[[i]] = T_star
    
    
    ## Confidence interval
    sigma_hat = f_Y^(-2)*H_hat
    alpha_sig = 0.05
    T_stat_abs = abs(unlist(T_stat, use.names=FALSE))
    q_star = quantile(T_stat_abs, 1-alpha_sig)
    CI = c(theta_hat-q_star*sigma_hat, theta_hat+q_star*sigma_hat)
    
    
    
    
    ## uncomment this for true theta setting
'    X = get_x(n, c)
    Y = get_y(n, theta, sig, 42)
    alpha = matrix(sort(rep(1/n, n)), nrow = n)
    theta_hat=theta
    H_hat = max(0.000001,abs((sum(alpha * (tau -   (Y <= theta_hat))))))
    e_multipliers = rnorm(n, 0, 1)
    T_star = -sum((H_hat^(-1/2) *((tau-(Y <= theta_hat)) * as.matrix(alpha)* e_multipliers)))
    T_stat[[i]] = T_star'
  }
  T_stats[[as.character(n)]]  = T_stat
  CIs[[as.character(n)]]  = CI
  d = density(unlist(T_stat), n=bootstraps)
  #d_norm = dnorm(unlist(T_stat))

  fn = glue(
    'qRF_location_',
    'q{formatC(tau*100, width=3, flag="0")}___',
    'new_var',
    'sigma{formatC(sig*100, width=3, flag="0")}___',
    'theta{formatC(theta*100, width=3, flag="0")}___',
    'n{formatC(as.integer(n), width=4, flag="0")}___',
    'bs{formatC(as.integer(bootstraps), width=4, flag="0")}',
    '.png'
  )
  
  png(file = fn)
  plot(d, ylab='', xlab='',ylim = c(0,0.45), main="Density plot of bootsrap test statistic")
  #ylim = c(0,1),
  x_std = rnorm(n, mean = 0, sd= 1)
  std_d = density(x_std, n=bootstraps)
  lines(std_d, col = 'blue')
  dev.off()
  
  #plot(seq(1,bootstraps,1), unlist(T_stats[["500"]]), ylab='', xlab='bootstrap', main="Test statistic vs standard normal points")
  #points(rnorm(100), col = 'blue')
  
}

 