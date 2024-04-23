#A function to count cluster sizes
counts = function(full, uniques)
{
  n = rep(0, length(uniques))
  for (i in seq(along.with = uniques))
  {
    n[i] = sum(full == uniques[i])
  }
  return(n)
}

#Create one hot encoded matrix C as in the paper
make_D <- function(beta_clust_labels, K) {
  p = length(beta_clust_labels)
  D <- matrix(0, nrow = p, ncol = K)
  for (l in 1:p) {
    for(k in 1:K){
      if(beta_clust_labels[l] == k)
      {
        D[l, k] <-1
      }
    }
  }
  return(D)
}
make_D_optim<- function(beta_clust_labels, K) {
  beta_clust_labels = as.numeric(beta_clust_labels)
  p <- length(beta_clust_labels)
  D <- matrix(0, nrow = p, ncol = K)
  for (l in 1:p) {
    D[l, beta_clust_labels[l]] <- 1
  }
  return(D)
}

make_D_dppm <- function(beta_clust_labels, K) {
  p = length(beta_clust_labels)
  D <- matrix(0, nrow = p, ncol = K)
  for (l in 1:p) {
    for(k in 1:K){
      if(beta_clust_labels[l] == k-1)
      {
        D[l, k] <-1
      }
    }
  }
  return(D)
}

#Calculate posterior probabilities of cluster indicators
get_clust_post <- function(prior_probs, log_lik) {
  post_probs_log <- log(prior_probs) + log_lik
  post_probs_log <- post_probs_log - max(post_probs_log)
  post_probs <- exp(post_probs_log) / sum(exp(post_probs_log))
  return(post_probs)
}

#Calculate number that appears maximum nunber of times in a vector
most_common_number<-function(x){
  # Get the counts of each unique element in the vector
  counts <- table(x)
  # Find the maximum count
  max_count <- max(counts)
  # Find the number that appears the maximum number of times
  most_common_number <- as.numeric(names(counts[counts == max_count]))
  return(most_common_number)
}



