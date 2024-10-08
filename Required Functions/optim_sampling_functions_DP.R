#sampling Functions DP

#Sample theta given clust_int, lambda, sigma_sqr, alpha
sample_theta_optim <- function(K, gamma_sqr, sigma_sqr, abun_dat_agg, resp, temp_beta_clust_labels, iter) {
  # Calculate G
  G <- (1 / gamma_sqr) * diag(K) + (1 / sigma_sqr) * crossprod(abun_dat_agg)
  
  # Perform Cholesky decomposition of G
  G_chol <- chol(G) #Gives the upper triangular matrix where G_chol'G_chol = G
  
  # Calculate theta_Sigma using the inverse of G
  theta_Sigma <- chol2inv(G_chol) # solve(G) or solve(G_chol)%*% solve(t(G_chol))
  
  #print(theta_Sigma )
  # Calculate theta_mean without matrix multiplication
  theta_mean = theta_Sigma %*% ((1/sigma_sqr)* crossprod(abun_dat_agg, as.matrix(resp))) 
  #print(theta_mean)
  # Generate theta_star by sampling from a multivariate normal distribution
  library(MASS)  # Make sure the MASS package is loaded
  theta_star <- cbind(mvrnorm(1, mu = theta_mean, Sigma = theta_Sigma))

  # Calculate H as a numeric table from sorted beta_clust_labels
  H <- rbind(as.numeric(table(sort(temp_beta_clust_labels))))
  #print(H)
  # q is always zero as per your comment
  q=0 # since 1^{T}\beta = 0
  # Calculate alp without matrix inversion
  #alp = solve(crossprod(t(G_chol)%*% t(H))) %*% (q - H %*% theta_star)
  alp = solve(H %*% theta_Sigma %*% t(H)) %*% (q - H %*% theta_star)
  
  # Calculate theta without matrix multiplication
  theta <- theta_star + theta_Sigma %*% t(H) %*% alp
  
  return(list(theta = theta, theta_star = theta_star))
}



#Sample sigma_sqr given clust_ind, sigma_sqr, alpha
sample_sigmasqr_optim<-function(a_sigma,n,b_sigma, resp, abun_dat_agg, theta)
{
  sigma_shape = (a_sigma + (n/2))
  #sigma_scale = b_sigma + 0.5* (t( resp- (abun_dat_agg %*% theta)) %*% ( resp- (abun_dat_agg %*% theta))) 
  sigma_scale = b_sigma + 0.5* crossprod(resp- (abun_dat_agg %*% theta))
  sigma_sqr = MCMCpack::rinvgamma(1, shape = sigma_shape, scale = sigma_scale)
  return(sigma_sqr)
}


#Sample gamma_sqr given theta, clust_ind, alpha, sigma_sqr
sample_gammasqr_optim<- function(a_gamma, K, theta)
{
  gamma_shape = a_gamma + (K-1)/2
  gamma_scale =  b_gamma + 0.5* crossprod(theta)
  gamma_sqr <-rinvgamma(1, shape = gamma_shape, scale = gamma_scale)
  return(gamma_sqr)
}

#Sample alpha
sample_alpha_optim<-function(alpha,p, a_alpha, b_alpha, K)
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

