rm(list = ls())
print(.libPaths())
libraries = c('grf', 'pROC', 'glue', 'caret', 'future.apply')
lapply(libraries,function(x)if(!(x %in% installed.packages())){
  install.packages(x)})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)
set.seed(42)


sink_on = function() sink(file = 'log.txt', append = TRUE)
sink_off = function() sink(file = NULL)

sink_on = function() return()
sink_off = function() return()
# Read in data ------------------------------------------------------------
df_org = read.csv('p2p.csv')
df_org$X = NULL

# Splitting into train and test data --------------------------------------
trainIndex = createDataPartition(df_org$status, p = 0.75, list = FALSE)

df_train = df_org[trainIndex, ]
df_test  = df_org[-trainIndex, ]

# Splitting data into target and covariates -------------------------------
Y = df_train$status
X = subset(df_train, select = -c(status))

Y_test = df_test$status
X_test = subset(df_test, select = -c(status))


# log-regression ----------------------------------------------------------
model <- glm(status ~ ., family=binomial(link='logit'), data = df_train)
summary(model)
# anova(model, test='Chisq')
Y_hat = predict(model)

sink_on()
print('Log-Regression')
print(pROC::roc(response = Y, predictor = Y_hat))
sink_off()

Y_hat = predict(model, newdata = subset(df_test, select=-c(status)), type = 'response')
roc = pROC::roc(response = Y_test, predictor = Y_hat)
plot(roc)

sink_on()
print('Test:')
print(roc)
sink_off()

# regression forest -------------------------------------------------------
sink_on()
print('')
print('GRF:')

regRF <- regression_forest(X = X, Y = Y, tune.parameters = 'all', num.trees = 5000)
print(regRF)
x = variable_importance(regRF, max.depth = 100)
row.names(x) = colnames(X)
print(round(x,4))
Y_hat = predict(regRF)

print("train:")
print(pROC::roc(response = Y, predictor = Y_hat$predictions))

Y_hat = predict(regRF, X_test)
print("test:")
print(pROC::roc(response = Y_test, predictor = Y_hat$predictions))
sink_off()
save(regRF, file = 'regression_forest.Rdata')
# get_Y_hat function ------------------------------------------------------

get_var_xs = function(X, var){
  return(seq(min(X[,var]), min(1, max(X[,var])), length.out = 100))
}
get_Y_hat_for_var = function(regRF = regRF, var = 'ratio001', i = NULL,
                             comp_variance = TRUE, X = X){
  x = get_var_xs(X = X, var = var)
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
var = 'ratio001'
Y_hat = get_Y_hat_for_var(
  regRF = regRF, var = var, i = 1,
  comp_variance = TRUE, X = X_test
)
save(Y_hat, file = 'ICE_one.Rdata')
png(file='ICE_one.png', bg = 'transparent')
  x = get_var_xs(X = X_test, var = var)
  sigma.hat = sqrt(Y_hat$variance.estimates)
  plot(x, Y_hat$predictions, ylim = range(Y_hat$predictions, 0, 1), 
       xlab = var, ylab = "ICE", type = "l", col = 'blue', lwd = 2)
  lines(x, Y_hat$predictions + 1.96 * sigma.hat, lty = 2, col ='blue')
  lines(x, Y_hat$predictions - 1.96 * sigma.hat, lty = 2, col ='blue')
dev.off()
# Mean ICE ----------------------------------------------------------------

plan(multisession) ## Run in parallel on local computer
nrow(X_test)
Y_hats = future_lapply(1:nrow(X_test), function(i) {
  get_Y_hat_for_var(
    regRF = regRF,
    var = 'ratio001',
    i = i,
    X = X_test,
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
  var = 'ratio001'
  x = get_var_xs(X = X_test, var = var)
  plot(x, Y_hat_predictions[,1], ylim = range(Y_hat_predictions, 0, 1), 
       xlab = var, ylab = "Mean ICE", type = "l", col = 'grey', lwd = 0.25)
  for(i in 2:min(ncol(Y_hat_predictions), 50000)){
    lines(x, Y_hat_predictions[, i], col='grey', lwd = 0.5)
  }
  lines(x, Y_hat$predictions, col='blue', lwd = 2)
  lines(x, Y_hat$predictions + 1.96 * sqrt(Y_hat$variance.estimates), lty = 2, col ='blue')
  lines(x, Y_hat$predictions - 1.96 * sqrt(Y_hat$variance.estimates), lty = 2, col ='blue')
dev.off()