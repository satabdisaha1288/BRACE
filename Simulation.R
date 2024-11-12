#Independent Data
library(rockchalk)
set.seed(333)
source("generate_compositional_data_replicates.R")
#Number of replicates
repl = 10
#Correlation coefficient
rho = 0.5
#Signal to noise ratio
snr = 1
# Number of predictors
p = 100
#Number of samples
n = 300
#Dependent structure 
beta = c(rep(-0.88,4), rep(-1.41,6), rep(-1.95,5), -1.16, 0.96,0, 0, 0, rep(1.04,6), 
         rep(0.51,4), rep(1.95,7))
#Add the zeroes
beta = c(beta, rep(0,p-length(beta)))
components = beta[which(beta!=0)]
# rho = rho[h]
sigma = mean(abs(components))/snr
simdata_dep<- generate_compositional_data_replicates(n = n, p = p, rho = rho,
                                                     sigma = sigma, beta = beta, m = 0.5 ,high.mag = 1:10,
                                                     beta0 = 1, intercept = FALSE, repl = repl, 
                                                     gamma = 2)
