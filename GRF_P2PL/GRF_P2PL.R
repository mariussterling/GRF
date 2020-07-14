rm(list = ls())
libraries = c('grf', 'pROC', 'glue', 'caret', 'future.apply')
lapply(libraries,function(x)if(!(x %in% installed.packages())){install.packages(x)})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)
set.seed(42)
log = ''

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

sink(file = 'log.txt')
print('Log-Regression')
print(pROC::roc(response = Y, predictor = Y_hat))
sink(file = NULL)

Y_hat = predict(model, newdata = subset(df_test, select=-c(status)), type = 'response')
roc = pROC::roc(response = Y_test, predictor = Y_hat)
plot(roc)

sink(file = 'log.txt')
print('Test:')
print(roc)
sink(file = NULL)

# regression forest -------------------------------------------------------
regRF <- regression_forest(X = X, Y = Y, tune.parameters = 'all')
Y_hat = predict(regRF)
sink(file = 'log.txt')
print("train:")
print(pROC::roc(response = Y, predictor = Y_hat$predictions))

Y_hat = predict(regRF, X_test)
print("test:")
print(pROC::roc(response = Y_test, predictor = Y_hat$predictions))
sink(file = NULL)
# plot --------------------------------------------------------------------

get_var_xs = function(X, var){
  return(seq(min(X[,var]), min(1, max(X[,var])), length.out = 100))
}
get_Y_hat_for_var = function(regRF = regRF, var = 'ratio001', i = NULL, comp_variance = TRUE, X = X){
  x = get_var_xs(X, var)
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

Y_hat = get_Y_hat_for_var(regRF = regRF, var = var, i = 1, X = X_test)
save(Y_hat, file = 'ICE_one.Rdata')
# ICE (Individual Conditional Expectation) --------------------------------
svg(file='ICE_one.svg', bg = "transparent")
  var = 'ratio001'
  x = get_var_xs(X, var)
  sigma.hat = sqrt(Y_hat$variance.estimates)
  plot(x, Y_hat$predictions, ylim = range(Y_hat$predictions, 0, 1), 
       xlab = var, ylab = "ICE", type = "l")
  lines(x, Y_hat$predictions + 1.96 * sigma.hat, col = 1, lty = 2)
  lines(x, Y_hat$predictions - 1.96 * sigma.hat, col = 1, lty = 2)
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
save(Y_hats, file = 'ICE_mean.Rdata')

Y_hat = data.frame(
  predictions = apply(sapply(Y_hats, function(x) x[, 1]), 1, mean),
  variance.estimates = apply(sapply(Y_hats, function(x) x[, 2]), 1, mean)
)