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
setwd("/BRACElet/Required Functions")
source("helperFunctions.R") # helper functions 
source("sampling.functions.DP.R") #sample theta, gamma_sqr, sigma_sqr, alpha
source("optim_sampling_functions_DP.R") # optimized version
#source("get_log_lik_dppm.R") # calculate log likehood for sampling cluster labels
source("get_log_lik_dppm_updated.R") # calculate log likehood for sampling cluster labels
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
        b_alpha = a_alpha / sqrt(p) 
        a_gamma = 0.001
        b_gamma = 0.001
        alpha_0 = 2 
        
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
          gamma_sqr[iter]<-sample_gammasqr_optim(a_gamma = a_gamma,b_gamma = b_gamma, K = n_clust, 
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
# type='structuredep0.5';
# setwd("/BRACElet/Output")
# datafile<-paste(type,'_snr',snr,'_n',n,'_p',p,'rho',rho,'repl',repl,"_outputs.RData",sep="")
# save(results,n,p,repl,rho, file=datafile)
