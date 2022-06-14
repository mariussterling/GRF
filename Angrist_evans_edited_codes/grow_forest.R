# This script was run with grf version 1.2.0.
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
rm(list = ls())
library('glue')
library(grf)
data = read.csv('angev80_recode_run1_line525.csv.xz')

FEATURES=data.frame(data$twoa.agem,
			  	data$twoa.agefstm,
			  	data$twoa.educm,
			  	as.numeric(data$twoa.blackm),
			  	as.numeric(data$twoa.hispm),
			  	as.numeric(data$twoa.othracem),
			  	twoa.incomed=round(data$twoa.incomed))
names(FEATURES)=1:ncol(FEATURES)

#labor income: twoa.incomem, worked for pay: twoa.workedm
DF.all=data.frame(
			  X=FEATURES,
			  Y=1 - as.numeric(data$twoa.workedm), # The outcome is whether the mother did not work
			  W=as.numeric(data$twoa.kidcount > 2),
			  I=as.numeric(data$twoa.samesex))

# remove ~15% of data with missing father's income
# roughly 4% of fathers have zero income after removing missing
#
# only consider married women
is.ok = !is.na(data$twoa.incomed) & (data$twoa.marital==0)
DF=DF.all[is.ok,]

#DF$Y.hat = predict(lm(Y ~ ., data = DF[,1:8]))


#W.lr = glm(W ~ ., data = DF[,c(1:7, 9)], family = binomial())
#DF$W.hat = predict(W.lr, type="response")



forest.iv = instrumental_forest(
  DF[,1:ncol(FEATURES)], DF$Y, DF$W, DF$I,
  min.node.size = 600,
  num.trees = 100000,
  ci.group.size = 125,
  sample.fraction = 0.2)

## epsilon_hat
#calculating E(Z) and E(WZ)

W.lr = regression_forest(DF[,c(1:7)], DF$W,  min.node.size = 800)

Y.lr = regression_forest(DF[,c(1:7)], DF$Y,   min.node.size = 800)

Z.lr = regression_forest(DF[,c(1:7)], DF$I,   min.node.size = 800)

DF$ZW = DF$I*DF$W
ZW.lr = regression_forest(DF[,c(1:7)], DF$ZW,   min.node.size = 800)

#save.image(file='myEnvironment.RData')

## only rerun from here after changing the test set

start_time <- Sys.time()

agefstm.vals= c(18,  22)
## qunatile grid
#incomed.vals = quantile(data$twoa.incomed, na.rm=TRUE, seq(0.025, 0.975, by = 0.05))

## uniform grid
up.lim = max(data$twoa.incomed, na.rm = TRUE)
lo.lim = min(data$twoa.incomed, na.rm = TRUE)
incomed.vals = seq(lo.lim, up.lim, length.out = 96)

dummy = expand.grid(AGEFSTM=agefstm.vals, INCOMED=incomed.vals)
X.test = data.frame(
  median(data$twoa.agem, na.rm=TRUE),
  dummy[,1],
  12,
  0, 0, 1,
  dummy[,2])
names(X.test)=1:ncol(X.test)

preds.iv = predict(forest.iv, X.test, estimate.variance = TRUE)
tau.hat = preds.iv$predictions
var.hat = preds.iv$variance.estimates
weights.hat = get_forest_weights(forest.iv, X.test)

names(X.test) <- paste0('X.', 1:(ncol(X.test)))

X.test$W.hat = predict(W.lr, newdata = X.test[,c(1:7)])$predictions
X.test$Y.hat = predict(Y.lr, newdata = X.test[,c(1:7)])$predictions
X.test$Z.hat = predict(Z.lr, newdata = X.test[,c(1:7)])$predictions
X.test$ZW.hat = predict(ZW.lr, newdata = X.test[,c(1:7)])$predictions



X.test$tau.hat = tau.hat

mu_hat = X.test$Y.hat - X.test$tau.hat*X.test$W.hat


## score function with two moments
psi_2 = sapply(1:nrow(X.test), function(k) (DF$Y-(X.test$tau.hat[k]*DF$W+mu_hat[k])))
psi_1 = sapply(1:nrow(X.test), function(k) DF$I*(DF$Y-(X.test$tau.hat[k]*DF$W+mu_hat[k])))

psi = list(psi_1, psi_2)
n = nrow(DF)


epsilon_hat=   lapply(1:nrow(X.test), function(j)
  -1* solve(matrix(c(X.test$ZW.hat[j],DF$W.hat[j], X.test$Z.hat[j], 1), nrow = 2, ncol = 2))
  %*% matrix(c(psi[[1]][,j], psi[[2]][,j]),nrow = 2, ncol =n))

## Taking only first element of the epsilon hat matrix
e_hat = sapply(epsilon_hat,'[',2,)


## test stat calculators
b = 500 #bootstraps
set.seed(100)
e_multipliers =  lapply(1:b, function(j) rnorm(n, 0, 1))
T_stat = sapply(1:b, function(j)
  ((var.hat^(-1/2)*weights.hat*t(e_hat))%*% e_multipliers[[j]])@x)
T_stat_abs =  abs(T_stat)
uniform_T_stat =  T_stat_abs
X.test$uniform_T = uniform_T_stat



alpha_sig = 0.05
## filtering the quantities by fixing mother's age 
X.test_18 = X.test[X.test$X.2 == 18, ]
X.test_22 = X.test[X.test$X.2 == 22, ]

uniform_T_max_18 =  apply(X.test_18$uniform_T, 1, max) # max of t_stats 
uniform_T_max_22 =  apply(X.test_22$uniform_T, 1, max) # max of t_stats 


uniform_q_star_18 =  quantile(uniform_T_max_18, 1-alpha_sig) # quantile of max_t_stat
uniform_q_star_22 =  quantile(uniform_T_max_22, 1-alpha_sig) # quantile of max_t_stat

output = data.frame(dummy, TAU=tau.hat, VAR=var.hat, 
                    Y = X.test$Y.hat,
                    W = X.test$W.hat,
                    Z = X.test$Z.hat,
                    ZW = X.test$ZW.hat,
                    T_stat = T_stat)

grids = nrow(X.test_18)
#write.table(output, file=glue("familysize_",grids,".out"))
write.table(output, file=glue("familysize_uniform",grids,".out"))
save(uniform_q_star_18, uniform_q_star_22, file =glue("critical_val_",grids,".Rdata"))

end_time <- Sys.time()


print(end_time - start_time)

