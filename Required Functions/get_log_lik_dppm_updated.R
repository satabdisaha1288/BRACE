get_log_lik_dppm_optim <- function(beta_clust_labels, p, j, K, X, sigma_sqr, gamma_sqr, y) {
  n <- length(y)
  log_lik <- numeric(K + 1)  # Initialize the vector

  # Calculate y' %*% y once since it's common to all cases
  yty <- crossprod(y)

  for (k in 0:K) {
    beta_clust_labels[j] <- k
    temp_beta_clust_labels <- beta_clust_labels[beta_clust_labels != 0]
    f = rbind(as.numeric(table(sort(as.numeric(temp_beta_clust_labels)))))
    n_clust <- length(table(temp_beta_clust_labels))

    if (length(temp_beta_clust_labels) == 0) {
      log_lik[k + 1] <- -0.5 * n * (log(2 * pi) + log(sigma_sqr)) -
        0.5 * (1/sigma_sqr) * yty
    } else {
      D <- make_D_optim(temp_beta_clust_labels, n_clust)
      X_sub <- X[, beta_clust_labels != 0]
      X_agg <- X_sub %*% D
      #Precision matrix
      G <- crossprod(X_agg) + (sigma_sqr / gamma_sqr) * diag(n_clust)
      # Partition A
      A_11 <- G[1:(n_clust-1), 1:(n_clust-1)]
      A_12 <- G[1:(n_clust-1), n_clust]
      a_n_clust <- G[n_clust, n_clust]

      # Define f_star
      f_star <- f[-n_clust]  # Example vector for f*
      f_K <- f[n_clust]

      # Calculate A_star
      A_star <- A_11 - ((1 / f_K) * (A_12 %*% t(f_star) + f_star %*% t(A_12))) + (a_n_clust / f_K^2) * (f_star %*% t(f_star))

      # Cholesky decomposition of G
      G_chol <- chol(A_star)
      # Calculate t(y) %*% X_agg efficiently
      #yt_X_agg <- crossprod(y, X_agg)
      # Calculate t(X_agg) %*% y efficiently
      # Calculate b

      X_agg_t_y <- crossprod(X_agg, y)
      b_star <- X_agg_t_y[1:(n_clust-1)]
      b_K <- X_agg_t_y[n_clust]
      b_tilde <- b_star - (b_K/f_K)*f_star
      # Calculate the determinant of G using Cholesky decomposition
      #log_det_G <- sum(log(diag(G_chol))) * 2
      # Calculate the intermediate term for log likelihood
      # Solve linear systems using Cholesky decomposition
      #G_inv_X_agg_t_y <- backsolve(G_chol, X_agg_t_y)
      G_inv_X_agg_t_y <- forwardsolve(t(G_chol), b_tilde)
      #G_inv_yt_X_agg <- backsolve(G_chol, t(X_agg) %*% y)

      # Calculate the intermediate term for log likelihood
      # intermediate_term <- (0.5 * sigma_sqr) * (
      #   crossprod(y) - sum(G_inv_yt_X_agg * G_inv_X_agg_t_y)
      # )
      intermediate_term <- (0.5 * (1/sigma_sqr)) * (
        crossprod(y) - crossprod(G_inv_X_agg_t_y)
      )
      # f_star <- f[-n_clust]  # Example vector for f*
      # f_K <- f[n_clust]
      # log_term4 <- 0.5 * log(sum(f_star^2) + f[n_clust]^2) - log(f[n_clust])
      # Calculate the log likelihood
      log_lik[k + 1] <-
        -(n / 2) * (log(2 * pi) + log(sigma_sqr))  +
        ((n_clust - 1) / 2) * (log(sigma_sqr) - log(gamma_sqr)) -
        (sum(log(diag(G_chol)))) - intermediate_term
    }
  }
  names(log_lik) <- 0:K
  return(log_lik)
}




