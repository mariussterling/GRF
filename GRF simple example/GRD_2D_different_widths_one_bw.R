## Checking coverage with different widths 

# Getting the path of your current open file ----
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

# Loading source files ----
source("support/packages_and_libraries.R") 
source("support/functions.R") 

ptm <- proc.time()
# Data initializing  ------------------------------------------------------
tau = c(0.5)
sig = 0.5
#widths = c(0.001,0.01, 0.04, 0.1, 0.2, 0.4)
widths = c(0.2)
c = 1
n1 = 100 #data points for x1 
b = 500 #number of bootstraps for MBS
reps = 20 #repititions of whole simulation
grids_x1 = 50 #grid points for CBs
set.seed(100)
node_size = 5 #bias controlling param of RFS
x2_fixed = c(0.3,0.5) #fix x2
alpha_sig = 0.05

rfs = list()
T_stats = list()
CIs = list()
CIs_std = list()
coverages = list()
coverages_std = list()
coverages_widths_PW = list()
coverages_widths_uniform = list()
T_stat = list()


# Data generation ----
set.seed(100)
for (width in widths){
## Sample
X = get_x(n1, seed=42) #no replications for X because X is deterministic
X = X[order(X$X1),] #ordering X wrt X1
theta_fun = function(X) theta_triangle(X, width)
theta_true = theta_triangle(X, width) + qnorm(tau)*sig
Y = get_y(X, theta_fun, sig, NULL,reps)
n = nrow(X)

## Test set
X_test = expand.grid(X1 = seq(-0.5,0.5,length.out=grids_x1), X2 = x2_fixed)
Y_test = get_y(X_test, theta_fun, sig, NULL,reps)
theta_true_test = theta_triangle(X_test, width) + qnorm(tau)*sig
grids = nrow(X_test)

## Fitting the forest ----
rf= fit_forest(X,Y,tau = tau,node_size = node_size)
theta_hat_test = predict_forest(rf,X_test)
w_test = weights_forest(rf,X_test)
mat <- do.call("cbind",theta_hat_test)
theta_hat_expected = rowMeans(mat)


## Calculations for test statistic ----
psi = lapply(1:reps, function(k) t(sapply(1:nrow(theta_hat_test[[k]]), function(j) tau - (Y[,k]<=theta_hat_test[[k]][j]))))
kde = lapply(1:reps, function (j) density(Y_test[,j], n=nrow(X))) #estimation of the density
f_Y= sapply(1:reps, function(j) unlist(approx(kde[[j]][["x"]], kde[[j]][["y"]], xout = c(theta_hat_test[[1]]))[2]))
V_hat = 1/f_Y
H_hat = sapply(1:reps, function(j) n* apply(w_test[[j]]*psi[[j]],1,var))
sigma_hat = (f_Y^(-2)*H_hat)^(1/2)
set.seed(5)
e_multipliers = lapply(1:reps, function(j) lapply(1:b, function(j) rnorm(n, 0, 1)))
T_stat = lapply(1:reps, function(k) sapply(1:b, function(j)
  (-(H_hat[,k]^(-1/2)*w_test[[k]]*psi[[k]])%*% e_multipliers[[k]][[j]])@x))

#T_stat = test_stat(X,Y,X_test,Y_test,theta_hat_test,w_test, tau= tau, b) 
T_stat_abs = lapply(1:reps, function(j) abs((T_stat[[j]])))

## just for simplicity, renaming test sets as original sets
# X = X_test
# Y = Y_test
# theta_hat = theta_hat_test
# theta_true = theta_true_test
# w = w_test

## Confidence interval with test stat ----
### alpha quantile of test statistic 
q_star = lapply(1:reps, function(j) rowQuantiles(T_stat_abs[[j]] , probs= c(1-alpha_sig)))
CI = confidence_interval(theta_hat_test, q_star, sigma_hat,reps)
## alpha quantile for standard normal
q_norm = lapply(1:reps, function(j) qnorm(1-alpha_sig))
CI_std = confidence_interval(theta_hat_test, q_norm, sigma_hat,reps)


## Calculating the coverage ----
### Coverage of true theta
coverage_true = coverage(theta_true_test, CI, grids, reps)
#coverage_std = coverage(theta_true_test, CI_std, grids, reps)
### Coverage of expected theta
#coverage_expected = coverage(theta_hat_expected, CI, grids, reps)
#coverage_expected_std = coverage(theta_hat_expected, CI_std, grids, reps)


## Calculation of uniform confidence bands ----
uniform_T_stat = lapply(1:reps, function(j) T_stat_abs[[j]]) 
uniform_T_max = lapply(1:reps, function(j) apply(uniform_T_stat[[j]], 1, max)) # max of t_stats 
uniform_q_star = lapply(1:reps, function(j) quantile(uniform_T_max[[j]], 1-alpha_sig)) # quantile of max_t_stat

uniform_CI = confidence_interval(theta_hat_test, uniform_q_star, sigma_hat, reps)
coverage_true_uniform = coverage(theta_true_test, uniform_CI, grids,reps)
#coverage_expected_uniform = coverage(theta_hat_expected, uniform_CI, grids,reps)

#
#----
coverages_widths_PW[[as.character(width)]] = coverage_true
coverages_widths_uniform[[as.character(width)]] = coverage_true_uniform
} 
# unlisting
coverages_widths_PW = as.numeric(unlist(coverages_widths_PW))
coverages_widths_uniform = as.numeric(unlist(coverages_widths_uniform))
# Plotting ----
## Plot of coverage with widths
plot(1)
png(file = glue('images/',
                'width_vs_coverage_',
                'n{formatC(as.integer(n),width=4, flag="0")}_',
                '.png'),
    width=1400, height=1400, res=300)
plot(widths,coverages_widths_PW,  
     type= 'l', xlab = bquote(eta), ylab = 'coverage',
     ylim=range( c(coverages_widths_PW, coverages_widths_uniform) ))
lines(widths,coverages_widths_uniform, col='magenta')
dev.off()
  
print(proc.time() - ptm)

# removing undesirable variables ----
rm(ptm)


