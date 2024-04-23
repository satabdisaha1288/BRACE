
#' Generate compositional data replicates
#'
#' @param n sample size
#' @param p number of compositional predictors
#' @param rho parameter used to generate the correlation matrix among predictors
#' @param sigma standard deviation for the noise terms, which are iid normal with mean 0
#' @param beta regression coefficients for the compositional predictors
#' @param m a scaler. For chosen components with large mean,
#'  log(p * m) is added to the non-normalized data 
#'  before the data are converted to compositional
#' @param high.mag an index vector with value(s) in [1,p], specifying which component(s) 
#' of compositions is has a mean of higher magnitude. 
#' @param beta0 coefficient of intercept
#' @param intercept logical value indicates whether to include an intercept
#' @param repl Number of replicates of the simulated data
#' @param gamma multiplicative constant during normalization
#'
#' @return A list simdata containing
#' @export
#'
#' @examples
generate_compositional_data_replicates<-function(n, p, rho,sigma, beta, m ,high.mag,
                                      beta0, intercept, repl, gamma)
{
  epsilon<-rnorm(n,0,sigma)
  theta <- rep(0, times = p)
  theta[high.mag] <- log(m * p)
  # Create the Sigma matrix using vectorized operations
  i <- 1:p
  j <- 1:p
  Sigma <- rho^abs(outer(i, j, "-"))
  #Generate matrix O
  data_mat <-replicate( repl, mvrnorm(n,theta, Sigma), simplify = FALSE)
  #  Transform and scale matrix 
  tran = lapply(data_mat, function(x) exp(gamma*x))
  rel.abun.dat = lapply(tran, function(x) x/rowSums(x))
  #Generate the compositional abundance matrix
  X<-lapply( rel.abun.dat, function(x) log(x))
  #Generate the responses
  y <- lapply( X, function(z) z %*% as.matrix(beta) + 
                 (as.integer(intercept) * beta0) + epsilon)
  simdata = list('y'=y, 'X_abs' = tran, intercept = intercept,
              'X_RelAbun' = rel.abun.dat,'X_logRelAbun'=X,'beta'=beta,'epsilon'=epsilon)
  return(simdata)
}
