# Example based on 
# https://github.com/grf-labs/grf/tree/master/experiments/

set.seed(1)
rm(list = ls())
libraries = c('grf', 'pROC', 'caret', 'future.apply', 'Hmisc',
              'lmtest', 'ggplot2', 'glue')
lapply(libraries,function(x)if(!(x %in% installed.packages())){
  install.packages(x)})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# Defining sink -----------------------------------------------------------
sink_on = function() sink(file = 'log.txt', append = TRUE)
sink_off = function() sink(file = NULL)

# sink_on = function() return()
# sink_off = function() return()


# Defining rsq and adj rsq ------------------------------------------------

rsq = function(x, x_pred){
  ss_tot = sum((x - mean(x))**2)
  ss_res = sum((x - x_pred)**2)
  return(1- ss_res/ss_tot)
}
rsq_adj = function(x, x_pred, p, model_with_const=TRUE){
  n = length(x)
  return(1 - (1-rsq(x, x_pred)) * (n - 1) / (n - p - 1 - ifelse(model_with_const, 1, 0)))
}

# Read in data ------------------------------------------------------------
data.all = read.csv("synthetic_data.csv")
data.all$schoolid = factor(data.all$schoolid)

DF = data.all[,-1]
school.id = as.numeric(data.all$schoolid)
school.mat = model.matrix(~ schoolid + 0, data = data.all)
school.size = colSums(school.mat)

# It appears that school ID does not affect pscore. So it can be ignored
# in modeling.
w.lm = glm(Z ~ ., data = data.all[,-3], family = binomial)
summary(w.lm)

# Setting up data set -----------------------------------------------------
trainIndex = createDataPartition(DF$Z, p = 0.75, list = FALSE)

W = DF$Z[trainIndex]
W_test = DF$Z[-trainIndex]
Y = DF$Y[trainIndex]
Y_test = DF$Y[-trainIndex]
X.raw = DF[,-(1:2)]

# Creation of dummy variables for multi label categorical variable C1 (student 
# race/ethnicity, 15 labels) and XC (School-lebel categorical variable for 
# urbanicity of the school)
C1_ = factor(X.raw$C1)
C1.exp = model.matrix(~ C1_ + 0)
XC_ = factor(X.raw$XC)
XC.exp = model.matrix(~ XC_ + 0)
X = cbind(X.raw[,-which(names(X.raw) %in% c("C1", "XC"))], C1.exp, XC.exp)
X_test = X[-trainIndex,]
X = X[trainIndex,]

# regression ----------------------------------------------------------
model <- lm(Y ~ ., data = cbind(X, W))
Y_hat = predict(model)

sink_on()
print('linear regression')
summary(model)
print(glue::glue("MSE: {mean((Y_hat - Y)**2)}"))
print(glue::glue("rsq: {round(rsq(Y,Y_hat),4)}"))
print(glue::glue("rsq_adj: {round(rsq_adj(Y,Y_hat, model$rank),4)}"))
sink_off()

Y_test_hat = predict(model, newdata = cbind(X_test, W=W_test))
sink_on()
print('Test:')
print(glue::glue("MSE: {mean((Y_test_hat - Y_test)**2)}"))
print(glue::glue("rsq: {round(rsq(Y_test, Y_test_hat),4)}"))
print(glue::glue("rsq: {round(rsq_adj(Y_test, Y_test_hat, model$rank),4)}"))
sink_off()

# regression forest -------------------------------------------------------
sink_on()
print('')
print('causal RF:')

regRF <- regression_forest(X = cbind(X, W), Y = Y, tune.parameters = 'all',
                         num.trees = 2000)
print(regRF)
x = variable_importance(regRF, max.depth = 50)
row.names(x) = colnames(X)
print(round(x,4))
Y_hat = predict(regRF)[ ,1]

print("train:")
print(glue::glue("MSE: {mean((Y_hat - Y)**2)}"))
print(glue::glue("rsq: {round(rsq(Y, Y_hat),4)}"))
print(glue::glue("rsq_adj: {round(rsq_adj(Y, Y_hat, ncol(regRF$X.orig)),4)}"))

Y_test_hat = predict(regRF, cbind(X_test, W=W_test))
print("test:")
print(glue::glue("MSE: {mean(((Y_test_hat - Y_test)**2)[,1])}"))
print(glue::glue("rsq: {round(rsq(Y_test, Y_test_hat),4)}"))
print(glue::glue("rsq_adj: {round(rsq_adj(Y_test, Y_test_hat, ncol(regRF$X.orig)),4)}"))
sink_off()
save(regRF, file = 'regression_forest.Rdata')

# get_Y_hat function ------------------------------------------------------
get_var_xs = function(X, var, l = 100){
  return(seq(min(X[,var]), max(X[,var]), length.out = l))
}
get_Y_hat_for_var = function(regRF = regRF, var, i = NULL,
                             comp_variance = TRUE, X = X, l =100){
  x = get_var_xs(X = X, var = var, l=l)
  X2 = data.frame(x)
  # X2 = data.frame(seq(-1, 1, length.out = 100))
  colnames(X2) = var
  if(is.null(i)){
    # X_fill = apply((X[,!colnames(X) %in% var]), 2, mean)
    X_fill = apply((X[,!colnames(X) %in% var]), 2, median)
  }else{
    X_fill = X[i, !colnames(X) %in% var]  
  }
  X2[, 2:(length(X_fill) + 1)] = X_fill
  colnames(X2)[2:ncol(X2)] = names(X_fill)
  
  Y_hat = predict(regRF, X2, estimate.variance = comp_variance)
  return(Y_hat)
}

# ICE (Individual Conditional Expectation) --------------------------------
var = 'S3'
Y_hat = get_Y_hat_for_var(
  regRF = regRF, var = var, i = 2,
  comp_variance = TRUE, X = cbind(X_test, W = W_test), l = 7
)
get_var_xs(X = cbind(X_test, W = W_test), var=var)
save(Y_hat, file = 'ICE_one.Rdata')
png(file='ICE_one.png', bg = 'transparent')
  x = get_var_xs(X = X_test, var = var, l = 7)
  sigma.hat = sqrt(Y_hat$variance.estimates)
  ci_l = Y_hat$predictions - 1.96 * sigma.hat
  ci_u = Y_hat$predictions + 1.96 * sigma.hat
  plot(x, Y_hat$predictions,
       xlab = var, ylab = "ICE", type = "l", col = 'blue', lwd = 2,
       # ylim = range(Y_hat$predictions, 0, 1),
       ylim = c(min(ci_l), max(ci_u)))
  lines(x, ci_u, lty = 2, col ='blue')
  lines(x, ci_l, lty = 2, col ='blue')
dev.off()
# Mean ICE ----------------------------------------------------------------

plan(multisession) ## Run in parallel on local computer
nrow(X_test)
Y_hats = future_lapply(1:nrow(X_test), function(i) {
  get_Y_hat_for_var(
    regRF = regRF,
    var = 'S3',
    i = i,
    X = cbind(X_test, W=W_test),
    l = 7,
    comp_variance = TRUE
  )},
  future.packages	= c('grf')
)
save(Y_hats, file = 'ICE_mean_test.Rdata')

# plan(multisession) ## Run in parallel on local computer
# Y_hats = future_lapply(1:nrow(X), function(i) {
#   get_Y_hat_for_var(
#     regRF = regRF,
#     var = 'ratio001',
#     i = i,
#     X = X,
#     comp_variance = TRUE
#   )},
#   future.packages	= c('grf')
# )
# save(Y_hats, file = 'ICE_mean_train.Rdata')

Y_hat = data.frame(
  predictions = apply(sapply(Y_hats, function(x) x[, 1]), 1, median),
  variance.estimates = apply(sapply(Y_hats, function(x) x[, 2]), 1, median)
)
Y_hat_predictions = sapply(Y_hats, function(x) x[,1])

png(file='ICE_mean.png', bg = "transparent")
var = 'S3'
  x = get_var_xs(X = X_test, var = var, l = 7)
  plot(x, Y_hat_predictions[,1], ylim = range(Y_hat_predictions, 0, 1), 
       xlab = var, ylab = "Mean ICE", type = "l", col = 'grey', lwd = 0.25)
  for(i in 2:min(ncol(Y_hat_predictions), 50000)){
    lines(x, Y_hat_predictions[, i], col='grey', lwd = 0.5)
  }
  lines(x, Y_hat$predictions, col='blue', lwd = 2)
  lines(x, Y_hat$predictions + 1.96 * sqrt(Y_hat$variance.estimates), lty = 2, col ='blue')
  lines(x, Y_hat$predictions - 1.96 * sqrt(Y_hat$variance.estimates), lty = 2, col ='blue')
dev.off()