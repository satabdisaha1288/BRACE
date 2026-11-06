#sampling Functions DP
#Sample theta given clust_int, lambda, sigma_sqr, alpha
#' Title Update value of theta
#'
#' @param K Number of clusters
#' @param gamma_sqr Variance of theta
#' @param sigma_sqr Variance of epsilon
#' @param abun_dat_agg Aggregated X matrix after multiplication with the cluster indicator matrix
#' @param resp Response vector y
#' @param temp_beta_clust_labels cluster labels for the regression coefficients ignoring the zero cluster
#' @param iter Iteration
#'
#' @return Updated value of theta
#' @return Updated value of compositional theta
#Sample theta given clust_int, lambda, sigma_sqr, alpha
sample_theta<-function(K,gamma_sqr, sigma_sqr, abun_dat_agg,resp, temp_beta_clust_labels)
{
  G = 1/(gamma_sqr)*diag(K) + (1/sigma_sqr)*t(abun_dat_agg) %*% abun_dat_agg
  #G_chol = chol(G)
  theta_Sigma = as.matrix(solve(G))
  theta_mean = theta_Sigma %*% ((1/sigma_sqr)* t(abun_dat_agg) %*% as.matrix(resp))
  #theta <- mvrnorm(1, mu= theta_mean, Sigma = theta_Sigma)
  theta_star <- cbind(mvrnorm(1, mu= theta_mean, Sigma = theta_Sigma))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Incorporate constrained sampling for theta using algorithm described in
  #Fast Simulation of Hyperplane-Truncated Multivariate Normal Distributions
  #muu = theta_mean
  #Sigma = theta_Sigma
  #x<-t(mvrnorm(1, mean = muu , Sigma = Sigma))
  #Get sorter cluster indices and form the constraint vector
  #H = cbind(1,1,1) # constraint coming from cluster assignment matrix
  #H = rbind(as.numeric(table(sort(beta_clust_labels[,iter][beta_clust_labels[,iter]!=0]))))
  H = rbind(as.numeric(table(sort(as.numeric(temp_beta_clust_labels)))))
  #H = t(as.matrix(rep(1, length(H)))) *(1/K)
  # constraint coming from cluster assignment vector beta_clust_labels[,iter]; updated at very iteration
  q = 0 # since 1^{T}\beta = 0
  alp = solve(H %*% theta_Sigma %*% t(H)) %*% (q - H %*% theta_star)
  theta = theta_star + theta_Sigma %*% t(H) %*% alp
  sum_beta = H%*%theta
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  return(list(theta = theta, theta_star = theta_star,
              sum_beta = sum_beta))
}

#Sample theta using elliptical slice sampling #Generic sampling

#Sample sigma_sqr given clust_ind, sigma_sqr, alpha
#Sample sigma_sqr given clust_ind, sigma_sqr, alpha
#' Title Update value of sigma_Sqr
#'
#' @param a_sigma Hyperpriors for the prior distribution of sigma
#' @param n Sample size
#' @param b_sigma Hyperpriors for the prior distribution of sigma
#' @param resp Response vector
#' @param abun_dat_agg Aggregated X matrix after multiplication with the cluster indicator matrix
#' @param theta Cluster means vector
#'
#' @return Updated value of sigmasqr
sample_sigmasqr<-function(a_sigma,n,b_sigma, resp, abun_dat_agg, theta)
{
  sigma_shape = (a_sigma + (n/2))
  sigma_scale = b_sigma + 0.5* (t( resp- (abun_dat_agg %*% theta)) %*% ( resp- (abun_dat_agg %*% theta)))
  sigma_sqr = MCMCpack::rinvgamma(1, shape = sigma_shape, scale = sigma_scale)
  return(sigma_sqr)
}


#Sample gamma_sqr given theta, clust_ind, alpha, sigma_sqr
#' Title Update value of gamma_sqr
#'
#' @param a_gamma Hyperprior for the prior distribution of gamma
#' @param b_gamma Hyperprior for the prior distribution of gamma
#' @param K Number of clusters
#' @param theta Cluster means vector
#'
#' @return Updated gamma_sqr value
sample_gammasqr<- function(a_gamma, b_gamma, K, theta)
{
  gamma_shape = a_gamma + (K-1)/2
  #gamma_shape = a_gamma + K/2 - 1
  gamma_scale = b_gamma + 0.5* (t(theta) %*% theta)
  gamma_sqr <-rinvgamma(1, shape = gamma_shape, scale = gamma_scale)
  return(gamma_sqr)
}


#Sample alpha
#' Update value of alpha
#'
#' @param alpha Previous iteration of alpha
#' @param p Number of predictors
#' @param a_alpha Hyperprior for the prior distribution of alpha
#' @param b_alpha Hyperprior for the prior distribution of alpha
#' @param K Number of clusters
#'
#' @return Updated value of alpha
sample_alpha<-function(alpha,p, a_alpha, b_alpha, K)
{
  #Auxillary variable
  omega<-rbeta(1, alpha + 1, p)
  a_star = (a_alpha + K - 1)
  b_star = p*(b_alpha - log(omega))
  pi_omega = a_star / (a_star + b_star)
  mix_ind <- rbinom(1, 1, pi_omega)
  if (mix_ind == 1) {
    this_shape <- a_alpha + K
  } else {
    this_shape <- a_alpha + K - 1
  }
  #this_scale = 1/ (b_alpha - log(omega))
  this_scale = (b_alpha - log(omega))
  alpha <- rgamma(1, shape = this_shape, scale = this_scale)
  return(alpha)
}

