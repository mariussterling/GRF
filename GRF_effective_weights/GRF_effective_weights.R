rm(list=ls())
library(grf)
library(ggplot2)
library(glue)
set.seed(42)


get_x = function(n){
  X = matrix((-(n / 2 - 1):(n / 2)) / (n / 2), nrow = n)
  colnames(X) = 'X1'
  return(X)
}
theta_triangle = function(x, width){
  pmax(1 - abs((x) / width), 0)
}

get_y = function(X, theta, sigma, seed=NULL){
  set.seed(NULL)
  n = nrow(X)
  return(theta(X) + rnorm(n, 0, sigma))
}

for (sig in c(0, 0.1)){
  rfs = list()
  for (n in c(500, 1000, 2000)){
    width = 0.2
    X = get_x(n)
    theta = function(X) theta_triangle(X, width)
    Y = get_y(X, theta, sig, 42)
    rf = regression_forest(X, Y)
    eps_tilde = Y - theta(X)
    alpha = get_sample_weights(rf)
    theta_tilde = theta(X) + alpha %*% eps_tilde
    
    pdf(file = glue('RF_theta_triangle___theta_tilde___sigma{formatC(sig*100, width=3, flag="0")}_n{formatC(n, width=4, flag="0")}.pdf'))
      plot(X[, 'X1'], Y, type = 'p', pch = 20, cex = 0.25, xlab='X1')
      lines(X[, 'X1'], theta_tilde, col = 'red')
      lines(X[, 'X1'], theta(X), col = 'blue', lty=2)
    dev.off()
    rfs[[as.character(n)]] = rf
  }
  
  x1s = get_x(100)
  alpha = lapply(rfs, function(x) get_sample_weights(x, newdata = x1s))
  rf_X.orig = lapply(rfs, function(x) x$X.orig[, 'X1'])
  for (i in 1:nrow(x1s)){
    pdf(file = glue('RF_theta_triangle___effective_weights___sigma{formatC(sig*100, width=3, flag="0")}___xi{formatC(x1s[i]*1000, width=4, flag="0")}.pdf'))
      # m = max(sapply(rfs, function(x)max(get_sample_weights(x)[])))
      plot(
        rf_X.orig[[1]],
        alpha[[1]][i, ],
        xlab = 'x', 
        ylab = parse(text = 'alpha[i](x)'),
        main = bquote('Effective weights'~alpha[i](x) ~'for' ~ X[i] == .(x1s[i])),
        type = 'l',
        ylim = c(0, 0.1)
      )
      lines(rf_X.orig[[2]], alpha[[2]][i, ], col = 'blue')
      lines(rf_X.orig[[3]], alpha[[3]][i, ], col = 'red')
    dev.off()
  }
}