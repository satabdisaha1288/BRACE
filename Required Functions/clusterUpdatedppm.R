
# functions to get the prior probabilities of 
# cluster membership
get_null_clust_prior <- function(beta_clust_labels, j, alpha_0, p) {
  #is_null <- as.numeric(beta_clust_labels == 0)
  is_null <- as.numeric(beta_clust_labels != 0)
  null_dict <- sort(table(is_null[-j]))
  null_keys <- as.vector(names(null_dict))
  
  # all 0 or all 1 need to adjust
  # All 0 or all 1 need to be adjusted
  if (identical(sort(setdiff(null_keys, c("0", "1"))), character(0))) {
    if ("0" %in% null_keys) {
      null_dict["1"] <- 0
    } else {
      null_dict["0"] <- 0
    }
    null_dict <- sort(null_dict)
  }
  #null_clust_values <- as.vector(null_dict)
  null_clust_values = as.vector(null_dict["0"])
  null_prior <- (null_clust_values + alpha_0 / 2) / (p - 1 + alpha_0)
  #print(paste0("null_prior =",null_prior))
  return(null_prior)
}

get_act_clust_prior <- function(beta_clust_labels, j, alpha) {
  is_null <- as.numeric(beta_clust_labels == 0)
  #temp_P <- sum(is_null[-j] == 0) + 1
  temp_P <- sum(is_null[-j] == 0) +  sum(is_null[j] == 0)
  clust_table_counts <- table(beta_clust_labels[-j][(beta_clust_labels[-j] != 0)])
  clust_table_counts <- sort(clust_table_counts)
  clust_table_counts <- as.vector(clust_table_counts)
  clust_prior <- clust_table_counts / (temp_P - 1 + alpha)
  clust_prior <- c(clust_prior, alpha / (temp_P - 1 + alpha))
  #print(paste0("clust_prior =",clust_prior))
  return(clust_prior)
}


get_clust_prior_dppm <- function(beta_clust_labels, alpha, alpha_0, p, j) {
  null_prior <- get_null_clust_prior(beta_clust_labels, j, alpha_0, p)
  clust_prior <- get_act_clust_prior(beta_clust_labels, j, alpha)
  # the second entry in the null-prior corresponds to the probability
  # of a variable being null
  prior_probs <- c(null_prior, clust_prior * (1-null_prior))
  return(prior_probs = prior_probs)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# function to take into account that everything can be null
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
get_prior_dppm <- function(j, beta_clust_labels, alpha_0, alpha, p) {
  # determining the null clusters and the clusters which are non-zero
  is_null <- as.numeric(beta_clust_labels == 0)
  temp_d_ind <- beta_clust_labels[(is_null == 0) & (1:length(beta_clust_labels) != j)]
  #temp_d_ind <- beta_clust_labels[(is_null[-j] == 0)]
  # it is possible for all predictors to be set to zero so we have
  # to deal with that in the prior probability calculation
  if (length(temp_d_ind) == 0) {
    # the prior if everything is null is the 2-component mixture 
    # model prior; the proposal is that the element will not be null
    K <- 1
    prior_probs <- get_null_clust_prior(beta_clust_labels, j, alpha_0, p)
    prior_probs <- c(prior_probs, 1-prior_probs)
    names(prior_probs)<-c(0:(K))
    return(list(prior_probs, K))
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # the case where there are non-null predictors
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # setup for prior probs
  clust_dict <- table(temp_d_ind)
  #clust_names <- as.vector(names(clust_dict))
  is_alone <- sum(beta_clust_labels[j] == temp_d_ind) == 0
  
  # rearrange the cluster indicators and lambdas
  K <- length(clust_names) + 1
  
  # get the prior probabilities
  prior_probs <- get_clust_prior_dppm(beta_clust_labels = beta_clust_labels, 
                                      alpha = alpha , alpha_0 = alpha_0, p = p, j = j)
  
  #names(prior_probs)<-c(0:(K))
  return(list(prior_probs,K))
  
}

update_beta_clust_labels_dppm <- function(beta_clust_labels, P_c, X, alpha, alpha_0,
                                          sigma_sqr, gamma_sqr,y) {
  
  #iter_seq <- sample(1:P_c, P_c, replace = FALSE)
  
  for (j in 1:P_c) {
    clust_names<-sort(unique(beta_clust_labels))
    K<-length(clust_names)
    table_counts = counts(beta_clust_labels, clust_names)
    names(table_counts)<-clust_names
    #which_table = which(customer[j] == tables)
    #which_clust_names = which(beta_clust_labels[j] == clust_names) # what is the jth coefficients cluster assignment?
    #which_clust_names = clust_names[beta_clust_labels[j] == clust_names] # what is the jth coefficients cluster assignment?
    which_clust_names = as.character(clust_names[beta_clust_labels[j] == clust_names])
    #table_counts[which_table] = table_counts[which_table] - 1
    table_counts[which_clust_names] = table_counts[which_clust_names] - 1 # remove the jth coefficient from the cluster
    #if (table_counts[which_table] == 0)
    #if the cluster becomes empty when the jth coefficient is removed, then remove that cluster
    if (table_counts[which_clust_names] == 0)
    {
      last_cluster = as.character(K-1)
      table_counts[which_clust_names] <- table_counts[last_cluster] # last cluster to replace this empty cluster
      #loc_z <- (beta_clust_labels == K) # who are in the last cluster
      loc_z <- (beta_clust_labels == K-1) # who are in the last cluster
      beta_clust_labels[loc_z] <- which_clust_names # move them up to fill the just emptied cluster
      table_counts <- table_counts[-K] # take out the last cluster, which is now empty
      #table_counts <- table_counts[-last_cluster] # take out the last cluster, which is now empty
      K = K - 1 # decrease the total number of clusters by 1
      #table_counts = table_counts[-which_clust_names]
      #clust_names = clust_names[-which_clust_names]
      #K = K - 1
    }
    #beta_clust_labels[j] <- -1 # ensures beta_clust_labels[j] does not get counted as a cluster
    #Get the prior probabilities
    result <- get_prior_dppm(j, beta_clust_labels, alpha_0, alpha, p)
    prior_probs <- result[[1]]
    #print(paste0("prior_probs = ", prior_probs))
    #K <- result[[2]]
    #print(paste0("K =" ,K))
    
    # get the log likelihood under each cluster
    log_lik <- get_log_lik_dppm_optim(beta_clust_labels, p ,j, K = K, X,sigma_sqr, gamma_sqr,y)
    #print(paste0("log_lik = ", log_lik))
    # normalize the posterior probs and sample new cluster membership
    post_probs <- get_clust_post(prior_probs, log_lik)
    #names(post_probs)<-1:(K+2)
    #print(paste0("post_probs = ", post_probs))
    # draw a sample of which cluster this new coefficient should belong to
    #new_clust_label <- sample(0:(K+1), size = 1,replace = TRUE, prob = post_probs)
    new_clust_label <- sample(0:K, size = 1,replace = TRUE, prob = post_probs)
    # rearranging again; need to account for if all is null
    is_null <- as.numeric(beta_clust_labels == 0) # gives index of non zero elements
    temp_d_ind <- beta_clust_labels[is_null == 0]
    
    #clust_names <- as.vector(names(table(temp_d_ind)))
    #K <- length(clust_names)
    
    # create a new cluster if necessary
    #if (new_clust_label == K+1)
    if (new_clust_label == K)
    {
      table_counts <- c(table_counts, 0)
      #K = K + 2
      K = K + 1
      names(table_counts)= c(0:(K-1))
    } 
    beta_clust_labels[j] <- new_clust_label
    new_clust_label = as.character(new_clust_label)
    #print(paste0("beta_clust_labels[j] =", beta_clust_labels[j]))
    table_counts[new_clust_label] <- table_counts[new_clust_label] + 1 # update the cluster at table counts
  }
  #print(paste0("beta_clust_labels =", beta_clust_labels))
  return(beta_clust_labels)
}

