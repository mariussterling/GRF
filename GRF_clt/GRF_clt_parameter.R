c_func = function(w, pi) {
  return(pi * log(1 - w) / log(w) / 2)
}
b_min = function(w, pi) {
  return(1 + log(1 - w) / (1 - log(w) / pi))
}


w = seq(0.01, 0.5, 0.01)
png(file='c_function.png', bg = "transparent")
  plot (w, c_func(w, 0.50), type = 'l', lwd = 2,
        ylab = expression('c('~omega~','~pi~')'),
        xlab = expression(omega))
  lines(w, c_func(w, 0.01), type = 'l', col = 2, lwd = 2)
  lines(w, c_func(w, 0.10), type = 'l', col = 3, lwd = 2)
  lines(w, c_func(w, 0.25), type = 'l', col = 4, lwd = 2)
dev.off()

png(file='beta_min_function.png', bg = "transparent")
  plot (w, b_min(w, 0.50), type = 'l', lwd = 2,
        ylab = expression(beta~min~'('~omega~','~pi~')'),
        xlab = expression(omega))
  lines(w, b_min(w, 0.01), type = 'l', col = 2, lwd = 2)
  lines(w, b_min(w, 0.10), type = 'l', col = 3, lwd = 2)
  lines(w, b_min(w, 0.25), type = 'l', col = 4, lwd = 2)
dev.off()
