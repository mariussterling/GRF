rm(list=ls())
library(grf)
library(ggplot2)
library(glue)
set.seed(42)

p = 1
n = 1000
sigma = 0.2

X = matrix(rnorm(n * p, 0.5, sigma), nrow = n)
colnames(X) = paste0('X', 1:p)

theta = function(x){
  x
}

get_y = function(X, sigma, seed=NULL){
  set.seed(NULL)
  n = nrow(X)
  return(theta(X) + rnorm(n, 0, sigma))
}


dens = list()
dens[['truth']] = density(theta(X))
for (sig in c(0.1, 0.2, 0.3)){
  Y = get_y(X, sig)
  rf = regression_forest(X, Y)
  eps_tilde = Y - theta(X)
  alpha = get_sample_weights(rf)
  theta_tilde = theta(X) + alpha %*% eps_tilde
  dens[[as.character(sig)]] = density(matrix(theta_tilde))
}

pdf('density_theta_tilde.pdf')
  plot(
    dens[['0.1']],
    col = 'black',
    xlab = glue(
      'X1\nN={n}, MeanBandwidth={round(mean(sapply(dens, function(x) x$bw)), 3)}'
    ),
    main = parse(text = 'Density~of~tilde(theta)(X)')
  )
  lines(dens[['0.2']], col = 'blue')
  lines(dens[['0.3']], col = 'red')
  lines(dens[['truth']], lty = 2)
dev.off()

png('density_theta_tilde.png')
  plot(
    dens[['0.1']],
    col = 'black',
    xlab = glue(
      'X1\nN={n}, MeanBandwidth={round(mean(sapply(dens, function(x) x$bw)), 3)}'
    ),
    main = parse(text = 'Density~of~tilde(theta)(X)')
  )
  lines(dens[['0.2']], col = 'blue')
  lines(dens[['0.3']], col = 'red')
  lines(dens[['truth']], lty = 2)
dev.off()