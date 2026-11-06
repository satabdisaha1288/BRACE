## ============================================================
## BRACE toy example: simulate compositional data, run the sampler,
## and compute the performance metrics reported in the paper.
##
## Run from the repository root:
##   source("toy_example.R")
## Packages: MASS, MCMCpack, salso
##
## With the default settings below this is a quick illustration: one
## simulated dataset (n = 300, p = 100, SNR = 1, AR(1) correlation) and a
## short chain of 300 iterations; runtime is about 3-4 minutes.
##
## To run the full simulation design of the paper, change the settings:
##   design = "block"           Scenario 1 (Section 3.7; Tables 1 and 3)
##   design = "ar1", rho = 0.5  AR(1) design (Supplementary Section S5.2)
##   p      in {100, 300, 1000}
##   snr    in {1, 5, 10}
##   n_rep  = 30
##   nburn  = 5000, nsamp = 3000   (8000 iterations, first 5000 discarded)
##
## Data: logistic-normal compositions; 34 nonzero coefficients in 8 clusters
## (columns 1-37) that sum to zero. 80/20 train/test split; the sampler is
## fit on the training set.
##
## Metrics (Sections 3.2-3.4 of the paper):
##   PE, L2 loss, TP / FP / FN, 95% predictive coverage and width,
##   ARI of the SALSO partition against the true clusters (nonzero
##   coefficients), number of SALSO clusters among the selected features,
##   and, when n_rep >= 2, the mean cross-replicate ARI.
## ============================================================

## ---- Settings ----------------------------------------------------------------
n      <- 300       # total sample size (80% used for training)
p      <- 100       # number of taxa
snr    <- 1         # signal-to-noise ratio
design <- "ar1"     # "ar1" or "block" (Scenario 1)
rho    <- 0.5       # AR(1) correlation, used when design = "ar1"
n_rep  <- 1         # number of simulated datasets (>= 2 gives cross-replicate ARI)
nburn  <- 150       # burn-in iterations
nsamp  <- 150       # post-burn-in iterations
seed   <- 333       # seed for data simulation
init_active <- 0.15 # share of predictors placed in an active cluster at the start

library(MASS)
library(MCMCpack)
library(salso)

fun_dir <- "Required Functions"
for (f in c("helperFunctions.R", "sampling.functions.DP.R", "optim_sampling_functions_DP.R",
            "get_log_lik_dppm_updated.R", "clusterUpdatedppm.R",
            "generate_compositional_data_replicates.R", "split_replicates_test_train.R"))
  source(file.path(fun_dir, f))

## ---- 1. Simulate data ------------------------------------------------------
set.seed(seed)
beta <- c(rep(-0.88, 4), rep(-1.41, 6), rep(-1.95, 5), -1.16, 0.96, 0, 0, 0,
          rep(1.04, 6), rep(0.51, 4), rep(1.95, 7))
beta  <- c(beta, rep(0, p - length(beta)))          # sum(beta) = 0
sigma <- mean(abs(beta[beta != 0])) / snr
truth <- beta != 0
true_clusters <- match(beta[truth], unique(beta[truth]))   # 8 true clusters

if (design == "ar1") {
  ## Sigma_ij = rho^|i - j|
  simdata <- generate_compositional_data_replicates(
    n = n, p = p, rho = rho, sigma = sigma, beta = beta, m = 0.5,
    high.mag = 1:10, beta0 = 1, intercept = FALSE, repl = n_rep, gamma = 2
  )
} else if (design == "block") {
  ## Scenario 1: within-cluster correlation 0.75 - 0.015|i - j|, between-cluster
  ## correlation 0.4 - 0.02|i - j| among nonzero coefficients, 0 otherwise,
  ## unit diagonal. Same transformation to compositions as the AR(1) design.
  d_ij  <- abs(outer(1:p, 1:p, "-"))
  same  <- outer(beta, beta, "==") & outer(truth, truth, "&")
  diffc <- outer(truth, truth, "&") & !same
  Sigma <- matrix(0, p, p)
  Sigma[same]  <- 0.75 - 0.015 * d_ij[same]
  Sigma[diffc] <- 0.4  - 0.02  * d_ij[diffc]
  diag(Sigma)  <- 1
  mu <- rep(0, p); mu[1:10] <- log(0.5 * p)
  X_rel <- replicate(n_rep, {
    Z <- exp(2 * mvrnorm(n, mu, Sigma))
    Z / rowSums(Z)
  }, simplify = FALSE)
  X_log <- lapply(X_rel, log)
  simdata <- list(X_RelAbun = X_rel, X_logRelAbun = X_log,
                  y = lapply(X_log, function(X) X %*% beta + rnorm(n, 0, sigma)))
} else stop("design must be \"ar1\" or \"block\"")

split <- split_replicates_test_train(split_ratio = 0.8,
                                     XlogRelAbun = simdata$X_logRelAbun,
                                     XRelAbun    = simdata$X_RelAbun,
                                     y           = simdata$y)

## ---- 2. Priors -----------------------------------------------------------------
a_sigma <- 0.001; b_sigma <- 0.001                  # sigma^2 ~ IG(a, b)
a_gamma <- 2.5;   b_gamma <- 1.5                    # gamma^2 ~ IG(2.5, 1.5)
a_alpha <- 1 / (0.75 * log(p))^2                    # alpha ~ Gamma(a, b)
b_alpha <- a_alpha / sqrt(p)
alpha_0 <- 2                                        # psi_0 ~ Beta(1, 1)
iterr   <- nburn + nsamp

## Adjusted Rand index
ari <- function(a, b) {
  tab <- table(a, b)
  comb2 <- function(x) sum(x * (x - 1) / 2)
  sum_ij <- comb2(tab); sum_a <- comb2(rowSums(tab)); sum_b <- comb2(colSums(tab))
  expected <- sum_a * sum_b / comb2(length(a))
  denom <- (sum_a + sum_b) / 2 - expected
  if (denom == 0) return(as.numeric(sum_ij == expected))   # both partitions trivial
  (sum_ij - expected) / denom
}

## ---- 3. Gibbs sampler for one replicate --------------------------------------
run_brace <- function(X, y) {
  n <- length(y)
  sigma_sqr <- gamma_sqr <- alpha <- numeric(iterr)
  sigma_sqr[1] <- 0.1; gamma_sqr[1] <- 0.1; alpha[1] <- 1
  beta_draws        <- matrix(0, p, iterr)
  beta_clust_labels <- matrix(NA, p, iterr)

  set.seed(1001)
  beta_clust_labels[, 1] <- sample(0:2, p, replace = TRUE,
                                   prob = c(1 - init_active, init_active / 2, init_active / 2))
  clust_names <<- 0:2                               # used by the sampler functions

  for (iter in 2:iterr) {
    if (iter %% 100 == 0) cat("Iteration", iter, "of", iterr, "\n")


    P_c <- p
    beta_clust_labels[, iter] <- update_beta_clust_labels_dppm(
      beta_clust_labels = beta_clust_labels[, iter - 1], P_c = P_c,
      alpha = alpha[iter - 1], X = X, alpha_0 = alpha_0, y = y,
      sigma_sqr = sigma_sqr[iter - 1], gamma_sqr = gamma_sqr[iter - 1]
    )

    labels  <- as.numeric(beta_clust_labels[, iter])
    active  <- labels != 0
    act_lab <- labels[active]
    K       <- length(table(act_lab))
    D       <- make_D_optim(act_lab, K)
    X_agg   <- X[, active, drop = FALSE] %*% D

    theta <- sample_theta(K = K, gamma_sqr = gamma_sqr[iter - 1],
                          sigma_sqr = sigma_sqr[iter - 1], abun_dat_agg = X_agg,
                          resp = y, temp_beta_clust_labels = act_lab)$theta
    sigma_sqr[iter] <- sample_sigmasqr_optim(a_sigma, n, b_sigma, y, X_agg, theta)
    gamma_sqr[iter] <- sample_gammasqr_optim(a_gamma, b_gamma, K, theta)
    alpha[iter]     <- sample_alpha_optim(alpha[iter - 1], p, a_alpha, b_alpha, K)

    beta_draws[active, iter] <- as.vector(D %*% theta)
  }

  post <- (nburn + 1):iterr
  list(B      = beta_draws[, post],
       labels = apply(beta_clust_labels[, post], 2, as.numeric),
       sigma_sqr = sigma_sqr[post])
}

## ---- 4. Metrics for one replicate ---------------------------------------------
brace_metrics <- function(fit, X_test, y_test) {
  B        <- fit$B
  PIP      <- rowMeans(fit$labels != 0)
  beta_hat <- rowMeans(B)
  selected <- PIP > 0.5 & rowMeans(abs(B) > 0.1) > 0.90     # selection rule, Section 3.3

  ## 95% posterior predictive intervals for the test set (Section 3.2)
  y_rep <- sapply(seq_len(ncol(B)), function(s)
    X_test %*% B[, s] + rnorm(length(y_test), 0, sqrt(fit$sigma_sqr[s])))
  pi_lo <- apply(y_rep, 1, quantile, 0.025)
  pi_hi <- apply(y_rep, 1, quantile, 0.975)

  ## SALSO point estimate of the partition, VI loss (Section 3.4)
  z_salso <- salso(t(fit$labels), loss = VI())

  list(
    PE       = mean((y_test - X_test %*% beta_hat)^2),
    L2       = sqrt(sum((beta_hat - beta)^2)),
    TP       = sum(selected & truth),
    FP       = sum(selected & !truth),
    FN       = sum(!selected & truth),
    Coverage = mean(y_test >= pi_lo & y_test <= pi_hi),
    Width    = mean(pi_hi - pi_lo),
    ARI      = ari(true_clusters, z_salso[truth]),
    Clusters = length(unique(z_salso[selected])),    # SALSO clusters of selected features
    sum_beta = max(abs(colSums(B))),
    z_salso  = z_salso,
    table    = data.frame(true_beta = beta, post_mean = beta_hat, PIP = PIP)
  )
}

## ---- 5. Run ------------------------------------------------------------------
res <- vector("list", n_rep)
for (r in seq_len(n_rep)) {
  cat(sprintf("\nReplicate %d of %d\n", r, n_rep))
  y      <- as.vector(split$y_split[[r]]$train_data)
  X      <- split$XlogRelAbun_split[[r]]$train_data
  y_test <- as.vector(split$y_split[[r]]$test_data)
  X_test <- split$XlogRelAbun_split[[r]]$test_data
  fit    <- run_brace(X, y)
  res[[r]] <- brace_metrics(fit, X_test, y_test)
}

## ---- 6. Report -----------------------------------------------------------------
metrics <- do.call(rbind, lapply(res, function(m)
  unlist(m[c("PE", "L2", "TP", "FP", "FN", "Coverage", "Width", "ARI", "Clusters")])))
cat("\nMetrics per replicate (true: 34 nonzero coefficients in 8 clusters)\n")
print(round(metrics, 3))
cat("Max |sum(beta)| over posterior draws:",
    signif(max(sapply(res, `[[`, "sum_beta")), 3), "\n")

if (n_rep >= 2) {
  pairs <- combn(n_rep, 2)
  cross_ari <- apply(pairs, 2, function(ij)
    ari(res[[ij[1]]]$z_salso[truth], res[[ij[2]]]$z_salso[truth]))
  cat(sprintf("Mean cross-replicate ARI (nonzero coefficients): %.3f\n", mean(cross_ari)))
}

print(round(res[[1]]$table[1:40, ], 2))
