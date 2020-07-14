rm(list = ls())
library(grf)
library(pROC)
library(glue)
set.seed(42)

df_org = read.csv('p2p.csv')
df_org$X = NULL


# Splitting into train and test data --------------------------------------
id_train = sort(sample(nrow(df_org), nrow(df_org)*.75))
#id_train = sort(sample(nrow(df_org), nrow(df_org)*.15))

df_train = df_org[id_train, ]
df_test  = df_org[setdiff(1:nrow(df_org), id_train), ]

Y = df_train$status
X = subset(df_train, select = -c(status, ratio037))
W = df_train$ratio037

Y_test = df_test$status
X_test = subset(df_test, select = -c(status, ratio037))
W_test = df_test$ratio037


# log-regression ----------------------------------------------------------
model <- glm(status ~ ., family=binomial(link='logit'), data = df_train)
summary(model)
anova(model, test='Chisq')
Y_hat = predict(model, newdata = subset(df_test, select=-c(status)), type = 'response')
#Y_hat = ifelse(Y_hat >= 0.1, 1, 0)

roc = pROC::roc(response = Y_test, predictor = Y_hat)
print(roc)
plot(roc)
print(roc$auc)



# regression forest -------------------------------------------------------
tau.forest <- regression_forest(X=cbind(X,W), Y=Y, tune.parameters = 'all')

Y_hat = predict(tau.forest)
roc = pROC::roc(response = Y, predictor = Y_hat$predictions)
print("train:")
print(roc)

Y_hat = predict(tau.forest, cbind(X_test,W_test))
roc = pROC::roc(response = Y_test, predictor = Y_hat$predictions)
print("test:")
print(roc)

# plot --------------------------------------------------------------------

var = 'ratio001'
X2 = data.frame(seq(min(X[,var]), min(1, max(X[,var])), length.out = 500))
# X2 = data.frame(seq(-1, 1, length.out = 100))
colnames(X2) = var
X2 = cbind(X2, cbind(X[, !grepl(var, colnames(X))], W)[2,])

tau_hat = predict(tau.forest, X2, estimate.variance = TRUE)
sigma.hat <- sqrt(tau_hat$variance.estimates)

plot(X2[, var], tau_hat$predictions, ylim = range(tau_hat$predictions, 0, 1), 
     xlab = var, ylab = "Y_hat", type = "l")
lines(X2[, var], tau_hat$predictions + 1.96 * sigma.hat, col = 1, lty = 2)
lines(X2[, var], tau_hat$predictions - 1.96 * sigma.hat, col = 1, lty = 2)


# Generalized Random forest -----------------------------------------------

tau.forest <- causal_forest(X=X, Y=Y, W=W, tune.parameters = 'all')

Y_hat = predict(tau.forest)
roc = pROC::roc(response = Y, predictor = Y_hat$predictions)
print("train:")
print(roc)

Y_hat = predict(tau.forest, X_test)
roc = pROC::roc(response = Y_test, predictor = Y_hat$predictions)
print("test:")
print(roc)


vi = grf::variable_importance(tau.forest, max.depth = 100)
row.names(vi) = colnames(X)
print(vi)


# plot --------------------------------------------------------------------

var = 'ratio001'
X2 = data.frame(seq(min(X[,var]), max(X[,var]), 0.01))
colnames(X2) = var
X2 = cbind(X2, X[, !grepl(var, colnames(X))][1,])
tau_hat = predict(tau.forest, X2, estimate.variance = TRUE)
sigma.hat <- sqrt(tau_hat$variance.estimates)

plot(X2[, var], tau_hat$predictions, ylim = range(tau_hat$predictions, -1, 2), 
     xlab = var, ylab = "CLATE", type = "l")
lines(X2[, var], tau_hat$predictions + 1.96 * sigma.hat, col = 1, lty = 2)
lines(X2[, var], tau_hat$predictions - 1.96 * sigma.hat, col = 1, lty = 2)

predict()
# causal forest -----------------------------------------------------------

tau.forest <- causal_forest(X, Y, W, tune.parameters = 'all')
Y_hat = predict(tau.forest, cbind(X_test, W_test))
roc = pROC::roc(response = Y_test, predictor = Y_hat$predictions)
print(roc)
plot(roc)
plot(Y_test, Y_hat$predictions, #ylim = range(Y_hat$predictions, 0, 2), 
     xlab = "x", ylab = "tau", type = "l")


# a -----------------------------------------------------------------------
tau.forest <- causal_forest(X, Y, W, num.trees = 4000)
tau.hat <- predict(tau.forest, X.test, estimate.variance = TRUE)
sigma.hat <- sqrt(tau.hat$variance.estimates)

plot(X.test[, 1], tau.hat$predictions, ylim = range(tau.hat$predictions + 1.96 * sigma.hat, tau.hat$predictions - 1.96 * sigma.hat, 0, 2), xlab = "x", ylab = "tau", type = "l")
lines(X.test[, 1], tau.hat$predictions + 1.96 * sigma.hat, col = 1, lty = 2)
lines(X.test[, 1], tau.hat$predictions - 1.96 * sigma.hat, col = 1, lty = 2)
lines(X.test[, 1], pmax(0, X.test[, 1]), col = 2, lty = 1)

# b -----------------------------------------------------------------------
tau.forest <- causal_forest(X, Y, W, tune.parameters = 'all')
# tau.forest.quantile = grf::quantile_forest(
#   X, Y, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975))


roc = apply(predict(tau.forest.quantile, X_test), 2, function(x){
  pROC::roc(response = Y_test, predictor = x)
})
plot(roc[[1]])
pROC::auc(roc)

roc = pROC::roc(response = Y_test, predictor = predict(tau.forest, X_test)[[1]])
plot(roc)
pROC::auc(roc)


# !!!VORLAGE --------------------------------------------------------------
#causal_forest(df)
# Generate data.
n <- 2000
p <- 10
X <- matrix(rnorm(n * p), n, p)
X.test <- matrix(0, 101, p)
X.test[, 1] <- seq(-2, 2, length.out = 101)

# Train a causal forest.
W <- rbinom(n, 1, 0.4 + 0.2 * (X[, 1] > 0))
Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
tau.forest <- causal_forest(X, Y, W)

# Estimate treatment effects for the training data using out-of-bag prediction.
tau.hat.oob <- predict(tau.forest)
hist(tau.hat.oob$predictions)

# Estimate treatment effects for the test sample.
tau.hat <- predict(tau.forest, X.test)
plot(X.test[, 1], tau.hat$predictions, ylim = range(tau.hat$predictions, 0, 2), xlab = "x", ylab = "tau", type = "l")
lines(X.test[, 1], pmax(0, X.test[, 1]), col = 2, lty = 2)

# Estimate the conditional average treatment effect on the full sample (CATE).
average_treatment_effect(tau.forest, target.sample = "all")

# Estimate the conditional average treatment effect on the treated sample (CATT).
# Here, we don't expect much difference between the CATE and the CATT, since
# treatment assignment was randomized.
average_treatment_effect(tau.forest, target.sample = "treated")

# Add confidence intervals for heterogeneous treatment effects; growing more trees is now recommended.
tau.forest <- causal_forest(X, Y, W, num.trees = 4000)
tau.hat <- predict(tau.forest, X.test, estimate.variance = TRUE)
sigma.hat <- sqrt(tau.hat$variance.estimates)
plot(X.test[, 1], tau.hat$predictions, ylim = range(tau.hat$predictions + 1.96 * sigma.hat, tau.hat$predictions - 1.96 * sigma.hat, 0, 2), xlab = "x", ylab = "tau", type = "l")
lines(X.test[, 1], tau.hat$predictions + 1.96 * sigma.hat, col = 1, lty = 2)
lines(X.test[, 1], tau.hat$predictions - 1.96 * sigma.hat, col = 1, lty = 2)
lines(X.test[, 1], pmax(0, X.test[, 1]), col = 2, lty = 1)

# In some examples, pre-fitting models for Y and W separately may
# be helpful (e.g., if different models use different covariates).
# In some applications, one may even want to get Y.hat and W.hat
# using a completely different method (e.g., boosting).

# Generate new data.
n <- 4000
p <- 20
X <- matrix(rnorm(n * p), n, p)
TAU <- 1 / (1 + exp(-X[, 3]))
W <- rbinom(n, 1, 1 / (1 + exp(-X[, 1] - X[, 2])))
Y <- pmax(X[, 2] + X[, 3], 0) + rowMeans(X[, 4:6]) / 2 + W * TAU + rnorm(n)

forest.W <- regression_forest(X, W, tune.parameters = "all")
W.hat <- predict(forest.W)$predictions

forest.Y <- regression_forest(X, Y, tune.parameters = "all")
Y.hat <- predict(forest.Y)$predictions

forest.Y.varimp <- variable_importance(forest.Y)

# Note: Forests may have a hard time when trained on very few variables
# (e.g., ncol(X) = 1, 2, or 3). We recommend not being too aggressive
# in selection.
selected.vars <- which(forest.Y.varimp / mean(forest.Y.varimp) > 0.2)

tau.forest <- causal_forest(X[, selected.vars], Y, W,
                            #W.hat = W.hat, Y.hat = Y.hat,
                            tune.parameters = "all")

# Check whether causal forest predictions are well calibrated.
test_calibration(tau.forest)
predict(tau.forest)$predictions
