get_log_lik_dppm_optim <- function(beta_clust_labels, p, j, K, X, sigma_sqr, gamma_sqr, y) {
  n <- length(y)
  log_lik <- numeric(K + 1)  # Initialize the vector
  
  # Calculate y' %*% y once since it's common to all cases
  yty <- crossprod(y)
  
  for (k in 0:K) {
    beta_clust_labels[j] <- k
    temp_beta_clust_labels <- beta_clust_labels[beta_clust_labels != 0]
    n_clust <- length(table(temp_beta_clust_labels))
    
    if (length(temp_beta_clust_labels) == 0) {
      log_lik[k + 1] <- -0.5 * n * (log(2 * pi) + log(gamma_sqr)) -
        0.5 * sigma_sqr * yty
    } else {
      D <- make_D_optim(temp_beta_clust_labels, n_clust)
      X_sub <- X[, beta_clust_labels != 0]
      X_agg <- X_sub %*% D
      #Precision matrix
      G <- crossprod(X_agg) + (sigma_sqr / gamma_sqr) * diag(n_clust)
      # Cholesky decomposition of G
      G_chol <- chol(G)
      # Calculate t(y) %*% X_agg efficiently
      #yt_X_agg <- crossprod(y, X_agg)
      # Calculate t(X_agg) %*% y efficiently
      X_agg_t_y <- crossprod(X_agg, y)
      # Calculate the determinant of G using Cholesky decomposition
      #log_det_G <- sum(log(diag(G_chol))) * 2
      # Calculate the intermediate term for log likelihood
      # Solve linear systems using Cholesky decomposition
      #G_inv_X_agg_t_y <- backsolve(G_chol, X_agg_t_y)
      G_inv_X_agg_t_y <- forwardsolve(t(G_chol), X_agg_t_y)
      #G_inv_yt_X_agg <- backsolve(G_chol, t(X_agg) %*% y)
      
      # Calculate the intermediate term for log likelihood
      # intermediate_term <- (0.5 * sigma_sqr) * (
      #   crossprod(y) - sum(G_inv_yt_X_agg * G_inv_X_agg_t_y)
      # )
      intermediate_term <- (0.5 * sigma_sqr) * (
        crossprod(y) - crossprod(G_inv_X_agg_t_y)
      )
      # Calculate the log likelihood
      log_lik[k + 1] <- 
        -(n / 2) * (log(2 * pi) + log(gamma_sqr)) +
        (n_clust / 2) * (log(sigma_sqr) - log(gamma_sqr)) -
        (sum(log(diag(G_chol)))) - intermediate_term
    }
  } 
  names(log_lik) <- 0:K
  return(log_lik)
}

get_log_lik_dppm <- function(beta_clust_labels, p ,j, K, X,sigma_sqr, gamma_sqr,y) 
{
  log_lik<- c()
  #for (k in 0:(K + 1)) 
  for (k in 0:(K)) {
    
    # set cluster indicator and find out how many are non-zero
    beta_clust_labels[j] <- k
    temp_beta_clust_labels <- beta_clust_labels[beta_clust_labels != 0]
    
    # if everyone is zero, log-lik is much simpler
    if (length(temp_beta_clust_labels) == 0) {
      log_lik[k+1] <- -(n/2)*(log(2*pi) + 
                                log(gamma_sqr)) -
        ((0.5*sigma_sqr)* (t(y)%*%y))} else {
          # the case where there are non-null predictors
          # getting the vectorized predictor matrix
          n_clust <- length(table(temp_beta_clust_labels))
          D <- make_D_optim(temp_beta_clust_labels, n_clust)
          X_sub <- X[, beta_clust_labels != 0]
          #X_agg <- (X_sub %*% D)[cc_ind, ]
          X_agg <- (X_sub %*% D)
          
          # precision matrix and cholesky
          G = t(X_agg) %*% X_agg + (sigma_sqr/gamma_sqr)*diag(n_clust)
          #g <- t(X_tilde) %*% ehat
          #G_chol <- chol(G)
          #gt_Ginv_g <- sum((backsolve(G_chol, g))^2)
          
          # log likelihood
          log_lik[k+1] <- -(n/2)*(log(2*pi) + log(gamma_sqr)) + 
            (n_clust/2)*(log(sigma_sqr) - log(gamma_sqr)) -
            (0.5 * log(det(G))) -((0.5*sigma_sqr)* (t(y)%*%y - 
                                                      (t(y) %*% X_agg) %*% solve(G) %*% (t(X_agg) %*% y)))
        }
  }
  names(log_lik)<-c(0:(K))
  return(log_lik)
}


