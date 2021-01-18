rm(list=ls())
library(grf)
library(ggplot2)
library(glue)
set.seed(42)

sig = 0
for (sig in c(0, 0.1)){
  # X = matrix(rnorm(n * p, 0.5, sigma), nrow = n)
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
  
  
  # n --------------------------------------------------------------------
  
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
    
    png(glue('RF_theta_triangle___theta_tilde___sigma{formatC(sig*100, width=3, flag="0")}_n{formatC(n, width=4, flag="0")}.png'))
      plot(X[, 'X1'], Y, type = 'p', pch = 20, cex = 0.25, xlab='X1')
      lines(X[, 'X1'], theta_tilde, col = 'red')
      lines(X[, 'X1'], theta(X), col = 'blue', lty=2)
    dev.off()
    rfs[[as.character(n)]] = rf
  }
  
  for (x1 in c(seq(0, 0.3, 0.025), seq(0.4, 1, 0.2))){
    png(glue('RF_theta_triangle___effective_weights___sigma{formatC(sig*100, width=3, flag="0")}___xi{formatC(x1*1000, width=4, flag="0")}.png'))
      m = max(sapply(rfs, function(x)max(get_sample_weights(x)[])))
      rf = rfs[[names(rfs)[1]]]
      i = which(min(abs(rf$X.orig - x1)) == abs(rf$X.orig - x1))[1]
      plot(
        rf$X.orig[-c(i), 'X1'],
        get_sample_weights(rf)[i, ][-c(i)],
        xlab = 'x', 
        ylab = parse(text = 'alpha[i](x)'),
        main = bquote('Effective weights'~alpha[i](x) ~'for' ~ X[i] == .(x1)),
        type = 'l',
        ylim = c(0, 0.1)
      )
      rf = rfs[[names(rfs)[2]]]
      i = which(min(abs(rf$X.orig - x1)) == abs(rf$X.orig - x1))[1]
      lines(rf$X.orig[-c(i), 'X1'], get_sample_weights(rf)[i, ][-c(i)], col = 'blue')
      rf = rfs[[names(rfs)[3]]]
      i = which(min(abs(rf$X.orig - x1)) == abs(rf$X.orig - x1))[1]
      lines(rf$X.orig[-c(i), 'X1'], get_sample_weights(rf)[i, ][-c(i)], col = 'red')
    dev.off()
  }
}