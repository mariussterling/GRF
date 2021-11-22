#################################
# This file contains all the libraries and packages required for running the
# simulations for constructing uniform confidence bands for GRF estimates

# Packages installation ----
#install.packages("matrixStats")
#install.packages("matrixStats")
#install.packages('pracma')
#install.packages('interp')

# Packages loading ----
rm(list=ls())
library(grf)
library(ggplot2)
library(glue)
library(rstudioapi)
library(parallel)
library(matrixStats)
set.seed(42)

rm(list=ls())
library(grf)
library(ggplot2)
library(glue)
library(rstudioapi)
library(parallel)
library(matrixStats)
library(interp)
set.seed(42)
