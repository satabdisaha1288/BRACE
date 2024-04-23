# Load the required libraries
library(MCMCpack)
library(stringr)
library(BhGLM)
library(salso)
#Load Input
simdata = simdata_dep
n = n
p = p
#Load required functions
setwd("/Users/satab/Documents/BRACElet/Required Functions")
source("helperFunctions.R") # helper functions 
source("sampling.functions.DP.R") #sample theta, gamma_sqr, sigma_sqr, alpha
source("optim_sampling_functions_DP.R") # optimized version
source("get_log_lik_dppm.R") # calculate log likehood for sampling cluster labels
source("clusterUpdatedppm.R") # sample cluster labels
source("split_replicates_test_train.R") # for splitting the data
train_test_split<-split_replicates_test_train(split_ratio = 0.8,
                                              XlogRelAbun = simdata$X_logRelAbun,                                             XRelAbun = simdata$X_RelAbun,
                                              y = simdata$y)
results<-list()
repp = 1
while(repp <= repl)
{
  skip_to_next <- FALSE
  tryCatch(
    {
      results<-list()
      n = length(train_test_split$y_split[[1]]$train_data)
      p = length(simdata$beta)
      repl = repl
      snr = 5
      rho = 0.5
      for(repp in 1: repl)
      {
        print(paste0("Starting replicate ",repp))
        y = train_test_split$y_split[[repp]]$train_data
        X = train_test_split$XlogRelAbun_split[[repp]]$train_data
        a_sigma = 0.001
        b_sigma = 0.001
        a_alpha = 1/(0.75*log(p))^2
        #a_alpha = 3#0.1 #3
        b_alpha = a_alpha / sqrt(p) #0.1
        a_gamma = 0.001
        b_gamma = 0.001
        alpha_0 = 2 #0.5
        
        #number of samples after burning
        nsamp = 5000
        # number of burn samples
        nburn = 5000
        iterr = nsamp + nburn
        #Initialize matrices and vectors
        theta<-list()
        abun_dat_agg<-list()
        sigma_sqr<-vector("numeric", length = iterr)
        sigma_sqr[1] <- 0.1
        gamma_sqr<-vector("numeric", length = iterr)
        gamma_sqr[1] <- 0.1
        alpha<-vector("numeric", length = iterr)
        alpha[1] = 1
        beta_clust_labels<- matrix(NA,nrow = p, ncol = iterr)
        #Set initial number of clusters
        #K = ifelse(round(alpha[1]*log(p)), round(alpha[1]*log(p)), 1)
        #K = 5
        #Dependent
        K = 3
        #Set clust names
        clust_names<-0:(K-1)
        #beta_clust_labels[,1] = c(rep(0,50), rep(1,5),rep(0,43),rep(2,2))
        set.seed(1001)
        beta_clust_labels[,1] = sample(0:2, p, prob = c(0.85, 0.10, 0.05), replace = TRUE)
        for (iter in 2:(nburn + nsamp)) 
        {
          print(paste0("Starting iteration ",iter))
          P_c = length(beta_clust_labels[,iter-1][which(beta_clust_labels[,iter-1]!=0)])
          beta_clust_labels[,iter] <-update_beta_clust_labels_dppm(beta_clust_labels = beta_clust_labels[, (iter - 1)], 
                                                                   P_c = P_c,alpha = alpha[iter - 1], X= X, alpha_0 = alpha_0,y=y,
                                                                   sigma_sqr = sigma_sqr[iter - 1], gamma_sqr = gamma_sqr[iter - 1])
          temp_beta_clust_labels <- beta_clust_labels[,iter][beta_clust_labels[,iter] != 0]
          #P_c = length(temp_beta_clust_labels)
          n_clust <- length(table(temp_beta_clust_labels))
          # make the cluster indicator matrix and vec covariate matrix
          D <- make_D_optim(beta_clust_labels = temp_beta_clust_labels, K = n_clust)
          X_sub <- X[, beta_clust_labels[,iter] != 0]
          abun_dat_agg[[iter]] <- (X_sub %*% D)
          theta[[iter]]<-sample_theta(K=n_clust, gamma_sqr = gamma_sqr[iter-1], sigma_sqr = sigma_sqr[iter - 1], 
                                      abun_dat_agg = abun_dat_agg[[iter]], resp=y, temp_beta_clust_labels = temp_beta_clust_labels)[[1]]
          sigma_sqr[iter]<-sample_sigmasqr_optim(a_sigma = a_sigma,n = n,b_sigma = b_sigma, 
                                                 resp = y, abun_dat_agg = abun_dat_agg[[iter]], theta = theta[[iter]])
          gamma_sqr[iter]<-sample_gammasqr_optim(a_gamma = a_gamma, K = n_clust, 
                                                 theta = theta[[iter]])
          alpha[iter]<-sample_alpha_optim(alpha = alpha[iter - 1], p = p, a_alpha = a_alpha, b_alpha = b_alpha, 
                                          K = n_clust)
        }
      }
    }, error = function(err) 
    { 
      print(paste("MY_ERROR:  ",err))
      skip_to_next <<- TRUE
    }
    
  )
  if(skip_to_next) { next }
  else 
  {
    # Return results for this replicate
    results[[repp]]<-list(
      beta_clust_labels = beta_clust_labels,
      abun_dat_agg =  abun_dat_agg,
      theta = theta,
      sigma_sqr = sigma_sqr,
      gamma_sqr = gamma_sqr,
      alpha = alpha
    )
    repp <- repp + 1
  }
}  

# #Dependent
# type='structuredep0.2';
# setwd("/Users/satab/Documents/BRACElet/Data")
# datafile<-paste(type,'_snr',snr,'_n',n,'_p',p,'rho',rho,'repl',repl,"_outputs.RData",sep="")
# save(results,n,p,repl,rho, file=datafile)


beta_clust_labels_train <- lapply(results, function(x) apply(x$beta_clust_labels, 2, as.numeric))

# Function to recreate a beta vector for a single iteration
recreate_beta_iteration <- function(iter, labels, theta) {
  zero_indices <- labels == 0
  beta_vector <- numeric(length(labels))
  non_zero_indices <- !zero_indices
  beta_vector[non_zero_indices] <- theta[labels[non_zero_indices]]
  return(beta_vector)
}

# Recreate new_beta_train using lapply and sapply
new_beta_train <- lapply(1:repl, function(repp) {
  sapply((nburn + 1):(iterr), function(iter) {
    recreate_beta_iteration(iter, beta_clust_labels_train[[repp]][, iter], results[[repp]]$theta[[iter]])
  })
})

new_beta_train<-lapply(new_beta_train, function(x) round(x,2))
#Create sum of the beta vector 
new_beta_sum<-sapply(new_beta_train, colSums)
#Plot histogram
hist(new_beta_sum, breaks = 20)
#Calculate probability of lying between (-0.1 to 0.1)
prob.zero.sum_new_beta<-lapply(sum_new_beta, function(x) 
  length(which(between(x[nburn+1: iterr], -0.2,0.2)==TRUE))/nsamp)
which(is.na(sum_new_beta))
# Mean of the sum of the beta vector for all replicates
mean_sum_new_beta = sapply(lapply(new_beta_train, colSums),mean)
# SD of the sum of the beta vector for all replicates
sd_sum_new_beta = sapply(lapply(new_beta_train, colSums),sd)

colmeans_list<-lapply(new_beta_train, rowMeans)
sd_new_beta = lapply(new_beta_train, function(x) apply(x,1,sd))
quantile_new_beta = lapply(new_beta_train, function(x) apply(x, 1, 
                                                             function(x) quantile(x, probs = c (0.025, 0.5,0.975))))
df<-list()

for(repp in 1:repl)
{
  df[[repp]]<-data.frame(beta = simdata$beta, beta_hat = quantile_new_beta[[repp]][2,], 
                         CI_low = quantile_new_beta[[repp]][1,],
                         CI_high =quantile_new_beta[[repp]][3,])
  df [[repp]]%>%
    insight::format_table(ci_brackets = c("(", ")")) %>%
    insight::export_table(format = "html")
}

#Calculate prediction error
PE<-list()
accuracy<-list()
a4<-list()
for(repp in 1: repl)
{
  yhat<-as.matrix(train_test_split$XlogRelAbun_split[[repp]]$test_data) %*% colmeans_list[[repp]]
  test_n = length(train_test_split$y_split[[repp]]$test_data)
  num1<-as.matrix(train_test_split$y_split[[repp]]$test_data) - yhat
  (PE[[repp]]<- (1/test_n) * crossprod(num1))
  accuracy[[repp]] = sqrt(sum((simdata$beta - as.matrix(colmeans_list[[repp]]))^2))
}
for(repp in 1:repl){
  y<- as.matrix(train_test_split$y_split[[repp]]$test_data)
  yhat<-as.matrix(train_test_split$XlogRelAbun_split[[repp]]$test_data) %*% colmeans_list[[repp]]
  a4[[repp]] <- measure.glm(y, yhat , family="gaussian")
}
a4.df<-do.call(rbind, a4)
a4.mean = colMeans(a4.df)
a4.sd = apply(a4.df, 2, sd)
#Prediction error
PE
#l2 loss
accuracy
(average_PE = mean(unlist(PE)))
(sd_PE = sd(unlist(PE)))
(average_accuracy = mean(unlist(accuracy)))
(sd_accuracy = sd(unlist(accuracy)))

RP = simdata$beta!=0
# PP4 = (samples4[,4]<0)|(samples4[,3]>0)
# FP4 = sum(PP4=="TRUE"&RP=="FALSE")
# TP4= sum(PP4=="TRUE"&RP=="TRUE")
# FN4 = sum(PP4=="FALSE"&RP=="TRUE")
# TN4 = sum(PP4=="FALSE"&RP=="FALSE")

#Calculate error rates
FP4 = TP4 = TN4 = FN4 <- vector()
for(repp in 1:repl)
{
  #which(RP == TRUE)
  a = ((round(quantile_new_beta[[repp]][1,],2) < -10^-6) & (round(quantile_new_beta[[repp]][3,],2) > 10^-6))
  b = (round(quantile_new_beta[[repp]][1,],2) == 0) & (round(quantile_new_beta[[repp]][3,],2) == 0)
  c = (round(quantile_new_beta[[repp]][2,],2) == 0)
  PP4 = a | b | c
  PP4 = !PP4
  which(PP4 == TRUE)
  FP4[repp] = sum(PP4=="TRUE"&RP=="FALSE")
  TP4[repp] = sum(PP4=="TRUE"&RP=="TRUE")
  FN4[repp] = sum(PP4=="FALSE"&RP=="TRUE")
  TN4[repp] = sum(PP4=="FALSE"&RP=="FALSE")
}
mean(FP4); sd(FP4)
mean(FN4); sd(FN4)


#Compute rand index
#true clustering
labels <- match(simdata$beta, unique(simdata$beta))
beta_salso <-lapply(beta_clust_labels_train, function(x) salso::salso(t(x), 
                                                                      maxNClusters = 10))

(arand.mean<-mean(sapply(beta_salso, function(x) pdfCluster::adj.rand.index(labels,x ))))
(arand.sd <- sd(sapply(beta_salso, function(x) pdfCluster::adj.rand.index(labels,x ))))



  