rm(list=ls())
library(grf)
library(dplyr)
library(ggplot2)
library(glue)
library(NNTbiomarker)
set.seed(42)


# Function definition for data --------------------------------------------
get_x = function(n, seed=NULL){
  set.seed(seed)
  X = data.frame(matrix(
    c(
      seq(-0.5, 0.5, length.out = n),
      round(runif(n = n, min = -0.5, max = 0.5), 3)
      ),
    nrow = n,
    ncol = 2
  ))
  colnames(X) = c('X1', 'X2')
  return(X)
}
get_x_grid = function(n){
  X = expand.grid(
    X1 = seq(-0.5, 0.5, length.out = n),
    X2 = seq(-0.5, 0.5, length.out = n)
  )
  return(X)
}

theta_triangle = function(x, width){
  pmax(1 - abs((x[, 'X1']) / width), 0)
}

get_y = function(X, theta, sigma, seed=NULL){
  set.seed(NULL)
  n = nrow(X)
  return(theta(X) + rnorm(n, 0, sigma))
}


# Computing ---------------------------------------------------------------
for (sig in c(0, 0.1)){
  rfs = list()
  for (n in c(500, 1000, 2000)){
    width = 0.2
    X = get_x(n, seed=42)
    theta = function(X) theta_triangle(X, width)
    Y = get_y(X, theta, sig, 42)
    rfs[[as.character(n)]] = rf = regression_forest(X, Y)
    eps_tilde = Y - theta(X)
    alpha = get_sample_weights(rf)
    theta_tilde = theta(X) + alpha %*% eps_tilde
    
    
    p = tibble(data.frame(cbind(X, Y))) %>% 
      ggplot(aes(X1, X2)) +
      geom_point(aes(colour = Y)) +
      scale_colour_gradient(low = "#132B43", high = "#56B1F7") +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      )
    fn = glue(
      'scatter/',
      'RF_theta_triangle___theta_tilde___',
      'sigma{formatC(sig*100, width=3, flag="0")}___',
      'n{formatC(n, width=4, flag="0")}_scatter',
      '.png'
    )
    ggsave(fn, plot = p, dpi = 300)
  }
  
  x1s = get_x_grid(101)
  alpha = lapply(rfs, function(x) get_sample_weights(x, newdata = x1s))
  alpha_max = max(sapply(alpha, max))
  rf_X.orig = lapply(rfs, function(x) x$X.orig)
  x1s_subset = expand.grid(c(0, 0.18, 0.2, 0.22, 0.9), c(0, 0.2, 0.9))
  for (n in as.integer(names(rfs))){
    for (i in 1:nrow(x1s_subset)){
      x1 = x1s_subset[i, 1]
      x2 = x1s_subset[i, 2]
      id = argmin(rowSums(cbind(rf_X.orig[[as.character(n)]][, 1] - x1, rf_X.orig[[as.character(n)]][, 2] - x2) ** 2))
      fn = glue(
        'raster/',
        'RF_theta_triangle___effective_weights___',
        'sigma{formatC(sig*100, width=3, flag="0")}___',
        'n{formatC(as.integer(n), width=4, flag="0")}___',
        'x_{formatC(x1*100, width=3, flag="0")}_{formatC(x2*100, width=3, flag="0")}',
        '.png'
      )
      p = tibble(cbind(x1s, alpha = alpha[[as.character(n)]][, id])) %>%
        ggplot(aes(X1, X2, fill = alpha)) +
        geom_raster(interpolate = FALSE) +
        scale_fill_gradient2(
          limits = c(0, alpha_max),
          name = c(expression(alpha[i](x[1], x[2])))) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
        ) +
        annotate(
          "point",
          x = rf_X.orig[[as.character(n)]][id, 1],
          y = rf_X.orig[[as.character(n)]][id, 2],
          pch = 4,
          colour = 'red',
          size = 3
        )
      
      ggsave(fn, plot = p, dpi = 300)
      fn = glue(
        'contour/',
        'RF_theta_triangle___effective_weights___',
        'sigma{formatC(sig*100, width=3, flag="0")}___',
        'n{formatC(as.integer(n), width=4, flag="0")}___',
        'contour___',
        'x_{formatC(x1*100, width=3, flag="0")}_{formatC(x2*100, width=3, flag="0")}',
        '.png'
       )
      
      breaks = c(seq(0, 0.07, 0.01), round(alpha_max, 2))
      contour = tibble(cbind(x1s, alpha = alpha[[as.character(n)]][, id])) %>%
        ggplot(aes(X1, X2, z = alpha)) +
        geom_contour_filled(
          aes(fill = stat(level)),
          breaks = breaks,
        ) +
        scale_fill_brewer(
          name = c(expression(alpha[i](x[1], x[2]))),
          palette = 'Blues',
          direction = 1,
          guide = 'legend',
          drop = FALSE
        ) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
        ) +
        annotate(
          "point",
          x = rf_X.orig[[as.character(n)]][id, 1],
          y = rf_X.orig[[as.character(n)]][id, 2],
          pch = 4,
          colour = 'red',
          size = 2
        )
      ggsave(fn, plot = contour, dpi = 300)
    }
  }
}

# https://ostechnix.com/create-video-pdf-files-linux/