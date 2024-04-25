#sampling Functions DP

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
sample_sigmasqr<-function(a_sigma,n,b_sigma, resp, abun_dat_agg, theta)
{
  sigma_shape = (a_sigma + (n/2))
  sigma_scale = b_sigma + 0.5* (t( resp- (abun_dat_agg %*% theta)) %*% ( resp- (abun_dat_agg %*% theta))) 
  sigma_sqr = MCMCpack::rinvgamma(1, shape = sigma_shape, scale = sigma_scale)
  return(sigma_sqr)
}


#Sample gamma_sqr given theta, clust_ind, alpha, sigma_sqr
sample_gammasqr<- function(a_gamma, K, theta)
{
  gamma_shape = a_gamma + K/2
  gamma_scale =  b_gamma + 0.5* (t(theta) %*% theta)
  gamma_sqr <-rinvgamma(1, shape = gamma_shape, scale = gamma_scale)
  return(gamma_sqr)
}


#Sample alpha
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

