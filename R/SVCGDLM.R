#' Temporal Correlation Matrix and Its Inverse
#'
#' Computes the temporal correlation matrix for a given length and decay parameter (phi), its inverse, and the log determinant.
#'
#' @param L Integer. The number of time points (must be positive).
#' @param phi Numeric. The temporal decay parameter (must be positive).
#'
#' @return A list with:
#'   \item{temporal_corr_inv}{Inverse of the temporal correlation matrix.}
#'   \item{log_deter}{Log determinant of the temporal correlation matrix.}
#'
#' @examples
#' temporal_corr_fun(5, 0.1)
#'
#' @export
#Temporal Correlation Function
temporal_corr_fun <- function(L, phi) {
  stopifnot(L > 0, phi > 0)  # Input validation
  L <- L + 1  # Adjust for 0-based indexing
  temporal_corr <- exp(-phi * abs(outer(1:L, 1:L, "-")))  # Vectorized
  temporal_corr_inv <- solve(temporal_corr)
  log_deter <- determinant(temporal_corr, logarithm = TRUE)$modulus
  return(list(temporal_corr_inv = temporal_corr_inv, log_deter = log_deter))
}

#' Update Polya-Gamma Latent Variables
#'
#' Generates random samples from the Polya-Gamma distribution for each observation in the Gibbs sampler.
#'
#' @param y Numeric vector. Outcome variable.
#' @param x Numeric matrix. Covariates for spatially varying coefficients.
#' @param z Numeric matrix. Covariates for fixed effects.
#' @param site_id Integer vector. Site or group IDs for each observation.
#' @param off_set Numeric vector (optional). Offset for each observation.
#' @param likelihood_indicator Integer. 0 for Bernoulli, 1 for Poisson/Negative Binomial.
#' @param log_r Numeric. Log of the r parameter (only for likelihood_indicator = 1).
#' @param theta Numeric vector. Coefficients for fixed effects.
#' @param Bs_old Numeric matrix. Current values of spatially varying coefficients.
#' @param eta_old Numeric vector. Current values for spatial random effects.
#'
#' @return A list with:
#'   \item{w}{Samples from the Polya-Gamma distribution.}
#'   \item{kse}{Pseudo-residuals for use in the Gibbs sampler.}
#'
#' @export
w_update <- function(y, x, z, site_id,off_set,likelihood_indicator, log_r, theta, Bs_old,eta_old) {
  n <- length(y)
  s <- nrow(Bs_old)
  off_set <- if (is.null(off_set)) rep(0, n) else off_set
  mean_w <- numeric(n)
  for (i in 1:s) {
    ids <- which(site_id == i)
    eta<-rep(eta_old[i],length(ids))
    mean_w[ids] <- off_set[ids] + 
      z[ids, ] %*% theta + 
      x[ids, ] %*% Bs_old[i, ]+eta
  }
  
  r <- exp(log_r)
  w <- numeric(n)
  kse <- numeric(n)
  
  if (likelihood_indicator == 0) {  # 0 for Bernoulli  and 1 for poisson
    w <- pgdraw(1, mean_w)
    kse <- (y - 0.50) / w
  }
  
  if (likelihood_indicator == 1) {
    w <- rpg.devroye(n,(y+r), mean_w)
    kse <- 0.50 * (y - r) / w
  }
  
  return(list(w = w, kse = kse))
}
#' Update Spatially Varying Coefficients
#'
#' Samples the spatially varying coefficients (Bs) in the Gibbs sampler.
#'
#' @param Bs_old Numeric matrix. Current values of Bs.
#' @param x Numeric matrix. Covariates for Bs.
#' @param z Numeric matrix. Covariates for fixed effects.
#' @param site_id Integer vector. Site/group ID for each observation.
#' @param off_set Numeric vector. Offset for each observation.
#' @param neighbors Numeric matrix. Adjacency matrix of site neighbors.
#' @param w Numeric vector. Polya-Gamma latent variables.
#' @param kse Numeric vector. Pseudo-residuals.
#' @param theta Numeric vector. Current fixed effect coefficients.
#' @param sigma2_beta Numeric. Prior variance for Bs.
#' @param gamma Numeric. Spatial correlation parameter.
#' @param eta Numeric vector. Spatial random effects.
#' @param corr_inv Matrix. Inverse of the temporal correlation matrix.
#'
#' @return A list with:
#'   \item{Bs}{Updated spatially varying coefficients.}
#'   \item{mu_beta}{Posterior mean of Bs.}
#'
#' @export
Bs_update <- function(Bs_old, x, z, site_id, off_set, neighbors, w, kse, theta, sigma2_beta, gamma, eta, corr_inv) {
  
  n <- length(w)  # sample size
  L <- ncol(x)-1    # number of columns in x
  s <- nrow(Bs_old)  # number of rows in beta_old
  w_mat <- matrix(w, n, (L+1))  # replicate w vector into matrix
  x_trans <- t(x)       # transpose of x
  Bs <- Bs_old    # initialize theta as beta_old
  mu_beta<-rmnorm (n=1,mean=rep(0, times=(L+1)),sigma2_beta*chol2inv(chol(corr_inv)))
  # Loop over each site
  for (i in 1:s) {
    mean_Bs_temp <- rep(0, (L+1))  # initialize mean_Bs_temp as a zero vector
    
    # Loop over each neighbor site
    for (k in 1:s) {
      mean_Bs_temp <- mean_Bs_temp + neighbors[i, k] * Bs[k, ]
    }
    
    ids <- which(site_id == i)  # find indices where site_id == j
    # Compute covariance matrix for Bs
    cov_Bs <- solve(
      x_trans[, ids] %*% (w_mat[ids, ] * x[ids, ]) +
        ((gamma * sum(neighbors[i, ]) + 1 - gamma) / sigma2_beta) * corr_inv)
    
    # Compute mean for Bs
    eta=rep(eta[i],length(ids))
    
    mean_Bs <- cov_Bs %*% (
      x_trans[, ids] %*% (w[ids]* (kse[ids]  -off_set[ids]- z[ids, ] %*% theta-eta)) +
        (corr_inv / sigma2_beta) %*% (gamma * mean_Bs_temp + (1 - gamma) * mu_beta)
    )
    
    # Sample Bs using multivariate normal distribution
    ind_norms <- matrix(rnorm(L+1), 1, (L+1))  
    Bs[i, ] <-as.matrix( mean_Bs + t(ind_norms %*% chol(cov_Bs)))
    
  }
  
  return(list(Bs=Bs,mu_beta=mu_beta))
}
#' Update Spatial Random Effects (Eta)
#'
#' Samples the spatial random effects for each site in the Gibbs sampler.
#'
#' @param Bs Numeric matrix. Current values of Bs.
#' @param x Numeric matrix. Covariates.
#' @param z Numeric matrix. Covariates for fixed effects.
#' @param site_id Integer vector. Site/group IDs.
#' @param off_set Numeric vector. Offset for each observation.
#' @param neighbors Numeric matrix. Adjacency matrix of site neighbors.
#' @param w Numeric vector. Polya-Gamma latent variables.
#' @param kse Numeric vector. Pseudo-residuals.
#' @param theta Numeric vector. Current fixed effect coefficients.
#' @param sigma2_eta Numeric. Prior variance for eta.
#' @param eta_old Numeric vector. Current values of eta.
#' @param rho Numeric. Spatial correlation parameter for eta.
#'
#' @return Numeric vector. Updated spatial random effects.
#'
#' @export
eta_update <- function(Bs, x, z, site_id, off_set, neighbors, w, kse, theta, sigma2_eta, eta_old, rho) {
  s <- length(eta_old)  # number of rows in eta_old
  eta <- eta_old  # initialize eta as eta_old
  
  for (i in 1:s) {
    ids <- which(site_id == i)  # find indices where site_id == j
    
    # Compute mean_eta_temp using neighbors
    mean_eta_temp <- sum(neighbors[i, ] * eta)
    
    # Compute variance 
    var_eta <- 1 / (sum(w[ids]) + (((rho * sum(neighbors[i, ]) + 1 - rho)) / sigma2_eta))
    
    # Compute mean 
    mean_eta <- var_eta * (
      sum(w[ids] * (kse[ids] - off_set[ids] - z[ids, ] %*% theta - x[ids, ] %*% Bs[i, ])) +
        (1 / sigma2_eta) * (rho * mean_eta_temp))
    
    # Sample eta[j] from a normal distribution
    eta[i] <- rnorm(1, mean = mean_eta, sd =sqrt (var_eta))
  }
  
  return(eta)
}
#' Gibbs Sampler for SVCGDLM
#'
#' Runs the Gibbs sampler for the Spatially Varying Coefficient Generalized Dynamic Linear Model.
#'
#' @param mcmc_samples Integer. Number of MCMC samples.
#' @param y Numeric vector. Outcome variable.
#' @param x Numeric matrix. Covariates for Bs.
#' @param z Numeric matrix. Covariates for fixed effects.
#' @param site_id Integer vector. Site/group IDs.
#' @param neighbors Numeric matrix. Adjacency matrix of site neighbors.
#' @param likelihood_indicator Integer. 0 for Bernoulli, 1 for Poisson/Negative Binomial.
#' @param off_set Numeric vector (optional). Offset.
#' @param theta Numeric vector. Initial fixed effect coefficients.
#' @param Bs_init Numeric matrix. Initial Bs.
#' @param eta_init Numeric vector. Initial eta.
#' @param sigma2_beta Numeric. Prior variance for Bs.
#' @param gamma Numeric. Spatial correlation parameter.
#' @param phi Numeric. Temporal correlation parameter.
#' @param sigma2_eta Numeric. Prior variance for eta.
#' @param rho Numeric. Spatial correlation parameter for eta.
#' @param log_r Numeric. Log of the r parameter (for likelihood_indicator = 1).
#'
#' @return A list with:
#'   \item{mcmc_samples}{Number of samples.}
#'   \item{Bs}{List of sampled Bs.}
#'   \item{mu_beta}{List of sampled mu_beta.}
#'   \item{eta}{List of sampled eta.}
#'   \item{w}{List of sampled Polya-Gamma variables.}
#'
#' @export
gibbs_sampler <- function(mcmc_samples,
                          y,
                          x,
                          z,
                          site_id,
                          neighbors,
                          likelihood_indicator,
                          off_set = NULL,
                          theta = NULL,
                          Bs_init=NULL,
                          eta_init=NULL,
                          sigma2_beta= NULL , 
                          gamma= NULL,
                          phi= NULL ,
                          sigma2_eta= NULL,
                          rho= NULL,
                          log_r = NULL) {
  
  L <- ncol(x)-1
  s <- ncol(neighbors)
  n <- length(y)
  # Initialize arrays for efficiency
  Bs <- mu_beta <- eta <- w <- vector("list", mcmc_samples)
  for (j in 1:mcmc_samples) {
    mu_beta_temp <- matrix(0, nrow = 1, ncol = (L+1))
    mu_beta[[j]] <- mu_beta_temp
  }
  for (j in 1:mcmc_samples) {
    Bs_temp <- matrix(0, nrow = s, ncol = (L+1))
    Bs[[j]] <- Bs_temp
  }
  for (j in 1:mcmc_samples) {
    eta_temp <- matrix(0, nrow = s, ncol = 1)
    eta[[j]] <- eta_temp
  }
  off_set <- if (is.null(off_set)) rep(0, n) else off_set
  Bs[[1]] <- if (is.null(Bs_init)) matrix(0, s, L + 1) else Bs_init
  eta[[1]] <- if (is.null(eta_init)) rep(0, s) else eta_init
  
  temporal_corr_info <- temporal_corr_fun(L, phi)
  D <- Matrix(diag(rowSums(neighbors)) - neighbors, sparse = TRUE)
  
  for (j in 2:mcmc_samples) {
    
    w_output <- w_update(y,
                         x,
                         z,
                         site_id,
                         off_set,
                         likelihood_indicator,
                         log_r,
                         theta,
                         Bs[[j - 1]],
                         eta[[j-1]])
    
    w[[j]] <- w_output[[1]]
    kse <- w_output[[2]]
    
    Bs_out <- Bs_update(Bs[[j - 1]],
                        x,
                        z,
                        site_id,
                        off_set,
                        neighbors,
                        w[[j]],
                        kse,
                        theta,
                        sigma2_beta,
                        gamma,
                        eta[[j - 1]],
                        temporal_corr_info[[1]])
    
    Bs[[j]]<- Bs_out$Bs
    mu_beta[[j]]<- Bs_out$mu_beta
    eta[[j]]=eta_update (Bs[[j]],
                         x,
                         z,
                         site_id,
                         off_set,
                         neighbors,
                         w[[j]], 
                         kse, 
                         theta,
                         sigma2_eta,
                         eta[[j-1]],
                         rho)
    
  }
  return(list(mcmc_samples=mcmc_samples,
              Bs = Bs,
              mu_beta=mu_beta,
              eta = eta,
              w=w))
  
}
#' Log Density of the Polya-Gamma Distribution
#'
#' Computes the log density of the Polya-Gamma distribution using a finite series expansion.
#'
#' @param x Numeric. Value(s) at which to evaluate the density.
#' @param b Numeric. Parameter for the Polya-Gamma distribution.
#' @param likelihood_indicator Integer. 0 for Bernoulli, 1 for Negative Binomial.
#'
#' @return Numeric. Log density value(s).
#'
#' @export
polya_gamma_log_density <- function(x, b = 1, likelihood_indicator = 0) {
  if (any(x <= 0)) return(ifelse(x <= 0, -Inf, NA))
  
  n_values <- 0:100  # Truncated series; consider dynamic adjustment
  pre_factor <- 1 / sqrt(2 * pi * x^3)
  log_pre_factor <- log(pre_factor)
  
  # Log of the constant term: (b - 1) * log(2) - log(Gamma(b))
  log_constant <- (b - 1) * log(2) - lgamma(b)
  
  # Compute log terms for the series
  log_terms <- lgamma(n_values + b) - lgamma(n_values + 1) + 
    log(2 * n_values + b) + log_pre_factor - 
    ((2 * n_values + b)^2) / (8 * x)
  
  # Numerically stable sum using log-sum-exp
  max_log_term <- max(log_terms)
  if (!is.finite(max_log_term)) return(-Inf)  # Handle overflow
  series_sum <- sum(exp(log_terms - max_log_term))
  if (series_sum <= 0) return(-Inf)
  
  log_series <- log(series_sum) + max_log_term
  log_density <- log_constant + log_series
  
  # Optional: Warn if series contribution is small
  if (log_series < -10) warning("Series sum may be underestimating for b = ", b, " and x = ", x)
  
  return(log_density)
}

#' Wrapper for Polya-Gamma Log Density
#'
#' Wrapper function to select correct log density for Bernoulli or Negative Binomial cases.
#'
#' @param x Numeric. Value(s) at which to evaluate the density.
#' @param y Numeric (optional). Outcome variable for Negative Binomial.
#' @param log_r Numeric (optional). Log r parameter for Negative Binomial.
#' @param likelihood_indicator Integer. 0 for Bernoulli, 1 for Negative Binomial.
#'
#' @return Numeric. Log density value(s).
#'
#' @export
polya_gamma_log_density_wrapper <- function(x, y = NULL, log_r = NULL, likelihood_indicator = 0) {
  if (likelihood_indicator == 0) {
    return(polya_gamma_log_density(x, b = 1))  # Bernoulli case, b = 1
  } else if (likelihood_indicator == 1) {
    if (is.null(y) || is.null(log_r)) stop("y and log_r must be provided when likelihood_indicator is 1")
    r <- exp(log_r)
    b <- y + r  # Negative binomial case
    return(polya_gamma_log_density(x, b = b))
  } else {
    stop("Invalid likelihood_indicator. Use 0 or 1.")
  }
}
#' Precompute Terms for Polya-Gamma Log Density
#'
#' Precomputes log factorials and powers for faster evaluation of the Polya-Gamma log density.
#'
#' @param n_terms Integer. Number of terms to precompute.
#'
#' @return A list of precomputed terms.
#'
#' @export
precompute_pg_log_terms <- function(n_terms = 200) {
  n <- 0:n_terms
  list(
    log_neg1_pow_n = log(abs((-1)^n)),  # Avoid complex numbers
    lgamma_n_plus_1 = lgamma(n + 1),     # Precompute log factorials
    n_values = n
  )
}

pg_log_terms <- precompute_pg_log_terms(200)  # Initialize once
#' Optimized Log Density for Polya-Gamma Distribution
#'
#' Computes the log density of the Polya-Gamma distribution using precomputed terms for speed.
#'
#' @param x Numeric. Value at which to evaluate the density.
#' @param b Numeric. Parameter for the Polya-Gamma distribution.
#' @param n_terms Integer. Number of terms in the series expansion.
#' @param precomp List. Precomputed terms (from \code{precompute_pg_log_terms}).
#'
#' @return Numeric. Log density value.
#'
#' @export
log_dpolya_gamma_optimized <- function(x, b, n_terms = 200, 
                                       precomp = NULL) {
  if(x <= 0) return(-Inf)
  if(b <= 0) stop("Parameter b must be positive")
  
  n <- 0:n_terms
  if (is.null(precomp)) {
    # Self-contained fallback
    log_terms <- log(abs((-1)^n)) + 
      lgamma(n + b) - 
      lgamma(n + 1) +
      log(2*n + b) - 
      0.5*log(2*pi*x^2)- 
      (2*n + b)^2/(8*x)
  } else {
    # Use precomputed terms
    log_terms <- precomp$log_neg1_pow_n + 
      lgamma(n + b) - 
      precomp$lgamma_n_plus_1 +
      log(2*n + b) - 
      0.5*log(2*pi*x^2)- 
      (2*n + b)^2/(8*x)
  }
  
  # Your superior sign handling
  max_log <- max(log_terms)
  signed_sum <- sum(sign((-1)^n) * exp(log_terms - max_log))
  
  if (signed_sum <= 0) return(-Inf)
  
  (b-1)*log(2) - lgamma(b) + max_log + log(signed_sum)
}

#' Vectorized Log Density for Polya-Gamma Distribution
#'
#' Computes the log density of the Polya-Gamma distribution for vectors of x and b.
#'
#' @param x Numeric vector. Values at which to evaluate the density.
#' @param b Numeric vector. Polya-Gamma parameters.
#' @param truncation Integer. Number of terms in the series expansion.
#'
#' @return Numeric vector. Log density values.
#'
#' @export
log_pg_density_vectorized <- function(x, b, truncation = 100) {
  # Vectorized input handling
  if (length(x) != length(b)) stop("x and b must have same length")
  
  # Initialize results
  results <- numeric(length(x))
  
  for (i in seq_along(x)) {
    if (x[i] <= 0) {
      results[i] <- -Inf
      next
    }
    if (b[i] <= 0) stop("All b values must be positive")
    
    # Rest of your calculation for single x[i], b[i] pair
    n <- 0:truncation
    log_const <- (b[i] - 1) * log(2) - lgamma(b[i])
    log_x_term <- -0.5 * (log(2*pi) + 3*log(x[i]))
    inv_8x <- 1/(8*x[i])
    
    log_terms <- lgamma(n + b[i]) - lgamma(n + 1) + 
      log(2*n + b[i]) + 
      log_x_term - 
      (2*n + b[i])^2 * inv_8x
    
    signs <- ifelse(n %% 2 == 0, 1, -1)
    max_log <- max(log_terms)
    signed_sum <- sum(signs * exp(log_terms - max_log))
    
    results[i] <- if (signed_sum <= 0) -Inf else 
      log_const + max_log + log(abs(signed_sum))
  }
  
  return(results)
}

#' Compute Q Function for MCEM Algorithm
#'
#' Computes the expected complete-data log-likelihood (Q function) in the MCEM algorithm for SVCGDLM.
#'
#' @param L Integer. Number of time points.
#' @param likelihood_indicator Integer. 0 for Bernoulli, 1 for Poisson/Negative Binomial.
#' @param param Numeric vector. Model parameters.
#' @param omega_mean Numeric vector. Mean Polya-Gamma variables.
#' @param Omega Matrix. Block-diagonal matrix of omega.
#' @param Bs_mean Numeric matrix. Mean of Bs samples.
#' @param mu_beta_mean Numeric vector. Mean of mu_beta samples.
#' @param eta_rep Numeric vector. Replicated eta.
#' @param chain_eta2 Numeric matrix. Eta samples.
#' @param chain_Bs2 Numeric array. Bs samples.
#' @param X_block Matrix. Block-diagonal design matrix.
#' @param site_id Integer vector. Site/group IDs.
#' @param y Numeric vector. Outcome variable.
#' @param z Numeric matrix. Covariates for fixed effects.
#' @param D Matrix. Difference matrix for spatial structure.
#'
#' @return Numeric. The value of the Q function.
#'
#' @export
Q_theta <- function(L=L,likelihood_indicator, param, omega_mean, Omega, Bs_mean, mu_beta_mean, eta_rep, chain_eta2, chain_Bs2, X_block, site_id, y, z, D) {
  # Extract parameters
  theta <- param[1:2]
  gamma <- param[3]
  rho <- param[4]
  phi <- param[5]
  sigma2_beta <- param[6]
  sigma2_eta <- param[7]
  
  n <- length(y)
  # Caching repeated calculations
  gamma_matrix <- Matrix(gamma * D + (1 - gamma) * diag(nrow(D)), sparse = TRUE)
  rho_matrix <- Matrix(rho * D + (1 - rho) * diag(nrow(D)), sparse = TRUE)
  inv_sigma_phi <- temporal_corr_fun(L, phi)[[1]]
  log_det <- temporal_corr_fun(L, phi)[[2]][[1]][[1]]
  
  # Precomputing determinants
  det_gamma <- determinant(gamma_matrix, logarithm = TRUE)$modulus
  det_rho <- determinant(rho_matrix, logarithm = TRUE)$modulus
  
  # Common calculations
  Q1 <- ((L + 1) / 2) * det_gamma - (s * (L + 1) / 2) * log(2 * pi) - (s * (L + 1) / 2) * log(sigma2_beta) - (s / 2) * log_det
  Q2k <- numeric(nrow(chain_Bs2))
  Q3 <- 0.5 * det_rho - (s / 2) * log(2 * pi) - (s / 2) * log(sigma2_eta)
  Q4k <- numeric(nrow(chain_Bs2))
  
  for (k in 1:nrow(chain_Bs2)) {
    Q4k[k] <- t(chain_eta2[, k]) %*% rho_matrix %*% chain_eta2[, k]
    dif <- as.vector(t(chain_Bs2[,, k]) - mu_beta_mean)
    Q2k[k] <- t(dif) %*% kronecker(gamma_matrix, inv_sigma_phi) %*% (dif)     
  }
  Q2 <- -1 / (2 * sigma2_beta) * mean(Q2k)
  Q4 <- -1 / (2 * sigma2_eta) * mean(Q4k)
  
  # Further calculations
  z_theta <- z %*% theta
  t_z_theta <- t(z_theta)
  XBS <- X_block %*% Bs_mean
  tXBS <- t(XBS)
  Omega_eta_rep <- Omega %*% eta_rep
  Omega_XBS <- Omega %*% XBS
  Omega_z_theta <- Omega %*% z_theta
  
  if (likelihood_indicator == 0) {
    pgdensity <- sum(sapply(omega_mean, function(omega) polya_gamma_log_density_wrapper(omega, likelihood_indicator = 0)))
    Q5 <- t(y - 0.5 - Omega %*% (XBS + eta_rep)) %*% z_theta - 0.5 * (t_z_theta %*% Omega_z_theta)
    Q6 <- -n * log(2) + t(y - 0.5) %*% (XBS + eta_rep)
    Q7 <- -0.5 * tXBS %*% Omega_XBS - 0.5 * t(eta_rep) %*% Omega_eta_rep - tXBS %*% Omega_eta_rep
    Q <- Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + pgdensity
  } else if (likelihood_indicator == 1) {
    r <- exp(param[8])
    pgdensity <- sum(log_pg_density_vectorized(omega_mean, y + r))
    
    Q5 <- sum(lgamma(y + r) - lgamma(r) - lgamma(y + 1)) - sum((y + r) * log(2))  
    Q6 <- 0.5 * t((y - r) - 2 * Omega %*% (XBS + eta_rep)) %*% z_theta - 0.5 * (t_z_theta %*% Omega_z_theta)
    Q7 <- 0.5 * t(y - r) %*% (XBS + eta_rep) - 0.5 * tXBS %*% Omega_XBS - 0.5 * t(eta_rep) %*% Omega_eta_rep - tXBS %*% Omega_eta_rep 
    Q <- Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + pgdensity
  } else {
    stop("Invalid likelihood_indicator value. Must be 0 or 1.")
  }
  
  return(as.numeric(Q))
}
#' Optimize SVCGDLM Parameters via MCEM
#'
#' Optimizes model parameters in the MCEM algorithm using parallel L-BFGS-B.
#'
#' @param likelihood_indicator Integer. 0 for Bernoulli, 1 for Poisson/Negative Binomial.
#' @param par0 Numeric vector. Initial parameter values.
#' @param nll Function. Negative log-likelihood function to be minimized.
#' @param cl Parallel cluster object.
#'
#' @return The output of \code{optimParallel}.
#'
#' @export
optimize_model <- function(likelihood_indicator, par0, nll, cl) {
  if (likelihood_indicator == 0) {
    par0 <- par0[1:7]  # Exclude r_init
    lower <- c(-Inf, -Inf, 0, 0, 0.0001, 0.001, 0.001)
    upper <- c(Inf, Inf, 0.999, 0.999, Inf, Inf, Inf)
  } else if (likelihood_indicator == 1) {
    par0 <- par0[1:8]  # Exclude r_init
    lower <- c(-Inf, -Inf, 0, 0, 0.0001, 0.001, 0.001, log(1))
    upper <- c(Inf, Inf, 0.999, 0.999, Inf, Inf, Inf, log(100))
  } else {
    stop("Invalid likelihood_indicator. Must be 0 or 1.")
  }
  
  # Optimize all parameters except r
  opt <- optimParallel(
    par = par0,  # Exclude r from continuous optimization
    fn = nll,
    lower = lower,
    upper = upper,
    method = "L-BFGS-B",
    control = list(maxit = 1000, factr = 1e7, pgtol = 1e-8, trace = FALSE),hessian=F,
    parallel = list(cl = cl)
  )
  
  return(opt)
}
#' Monte Carlo EM Algorithm for SVCGDLM
#'
#' Runs the Monte Carlo Expectation-Maximization algorithm for the SVCGDLM.
#'
#' @param mcmc_samples Integer. Number of MCMC samples.
#' @param y Numeric vector. Outcome variable.
#' @param x Numeric matrix. Covariates for Bs.
#' @param z Numeric matrix. Covariates for fixed effects.
#' @param site_id Integer vector. Site/group IDs.
#' @param neighbors Numeric matrix. Adjacency matrix of site neighbors.
#' @param likelihood_indicator Integer. 0 for Bernoulli, 1 for Poisson/Negative Binomial.
#' @param off_set Numeric vector (optional). Offset.
#' @param theta_init Numeric vector. Initial fixed effect coefficients.
#' @param Bs_init Numeric matrix. Initial Bs.
#' @param eta_init Numeric vector. Initial eta.
#' @param sigma2_beta_init Numeric. Initial prior variance for Bs.
#' @param gamma_init Numeric. Initial spatial correlation parameter.
#' @param phi_init Numeric. Initial temporal correlation parameter.
#' @param sigma2_eta_init Numeric. Initial prior variance for eta.
#' @param rho_init Numeric. Initial spatial correlation parameter for eta.
#' @param r_init Numeric. Initial value for r (Negative Binomial).
#' @param tol Numeric. Convergence tolerance.
#' @param max_iter Integer. Maximum number of EM iterations.
#' @return A list with parameter estimates, standard errors, predictions, and MCMC samples.
#' @export 
MCEM <- function(mcmc_samples,
                 y,
                 x,
                 z,
                 site_id,
                 neighbors,
                 likelihood_indicator,
                 off_set = NULL,
                 theta_init = NULL, 
                 Bs_init = NULL,
                 eta_init = NULL,
                 sigma2_beta_init= NULL , 
                 gamma_init= NULL,
                 phi_init= NULL ,
                 sigma2_eta_init= NULL,
                 rho_init= NULL,
                 r_init = NULL,
                 tol,
                 max_iter) {
  
  # create x block
  L <- ncol(x)-1
  n <- length(y)
  xblock<-lapply(1:s,function(i) x[site_id==i,])
  X_block <- bdiag(xblock)
  #MCEM algorithm
  par0_old=c(theta_init,gamma_init,rho_init,phi_init,sigma2_beta_init,sigma2_eta_init,r_init)
  par0_new=par0_old
  Bs_new=NULL
  eta_new=NULL
  for (iter in 1:max_iter) {
    #E-step: Gibbs sampling
    gibbs_result <-gibbs_sampler(mcmc_samples,
                                 y,
                                 x,
                                 z,
                                 site_id,
                                 neighbors,
                                 likelihood_indicator,
                                 off_set,
                                 theta=par0_new[1:2], 
                                 Bs_init = Bs_new,
                                 eta_init = eta_new,
                                 sigma2_beta=par0_new[6], 
                                 gamma=par0_new[3],
                                 phi=par0_new[5],
                                 sigma2_eta=par0_new[7],
                                 rho=par0_new[4],
                                 log_r=log(par0_new[8]))   
    
    m<-mcmc_samples
    burn_in<-0.2*m
    thin<-5
    #BS and mu_beta chain
    chain_Bs<-simplify2array(gibbs_result[[2]])
    chain_Bs2<-chain_Bs[,,seq(burn_in+1,m,thin)]
    Bs_means<-apply(chain_Bs2, c(1, 2), mean) 
    Bs_sd<-apply(chain_Bs2, c(1, 2), sd)
    # Calculate 0.025 and 0.975 quantiles
    Bs_quantile_025 <- apply(chain_Bs2, c(1, 2), quantile, probs = 0.025)
    Bs_quantile_975 <- apply(chain_Bs2, c(1, 2), quantile, probs = 0.975)
    Bs_mean<-as.vector(Bs_means)
    Bs_new<-Bs_means
    #mu_beta
    chain_mu_beta<-simplify2array(gibbs_result[[3]])
    chain_mu_beta2<-chain_mu_beta[,seq(burn_in+1,m,thin)]
    mu_beta_mean<-rowMeans(chain_mu_beta2)
    #eta chain
    chain_eta<-simplify2array(gibbs_result[[4]])
    chain_eta2<-chain_eta[,seq(burn_in+1,m,thin)]
    eta_mean<-rowMeans(chain_eta2)
    eta_rep=unlist(lapply(1:s,function(i) rep(eta_mean[i],sum(site_id==i))))
    eta_new=eta_mean
    #omega chain
    omega_matrix <- do.call(rbind, gibbs_result[[5]][seq(burn_in+1,m,thin)])
    omega_mean <- colMeans(omega_matrix)
    Omega <- bdiag(lapply(split(omega_mean, site_id), diag))
    Omega <- as.matrix(Omega)
    
    #M-step: Update fixed effects and random effect 
    nll <- function(x) -1 * Q_theta(L=L,likelihood_indicator=likelihood_indicator,param = x,omega_mean=omega_mean, Omega = Omega,Bs_mean=Bs_mean,mu_beta_mean=mu_beta_mean,eta_rep=eta_rep,
                                    chain_eta2=chain_eta2,chain_Bs2=chain_Bs2,X_block=X_block, site_id = site_id, y = y, 
                                    z = z, D = D)
    
    # Optimize using L-BFGS-B with optimParallel
    opt <- optimize_model(likelihood_indicator = likelihood_indicator, par0 = par0_new, nll = nll, cl = cl)
    # update parameters
    if (length(opt$par)==7)  {
      par0_new=round(opt$par, 3)
      covariance_matrix <- solve(opt$hessian)
      std_errors <- sqrt(diag(covariance_matrix))
    }
    if(length(opt$par)==8) {
      par0_new=round(c(opt$par[1:7],exp(opt$par[8])), 3)
    }
    cat("par0_new:",par0_new,"\n")
    if (max(abs(par0_old-par0_new)) < tol) {
      
      par0_old=par0_new
      
      break
    }
    par0_old=par0_new
    
  }
  
  if (likelihood_indicator == 0) {
    logit_p<-rep(0, times=n)
    for(j in 1:s){
      logit_p[site_id == j]<-z[site_id == j,]%*%par0_old[1:2] + 
        x[site_id == j,]%*%Bs_new[j,]+eta_new[j]
    }
    
    phat<-exp(logit_p)/(1 + exp(logit_p))
    
    return(list(par0_old=par0_old,
                std_errors=std_errors,
                prediction=phat,
                Bs_means=Bs_means,
                chain_Bs=chain_Bs,
                Bs_sd=Bs_sd,
                Bs_quantile_025=Bs_quantile_025,
                Bs_quantile_975=Bs_quantile_975
    ))
  }
  if (likelihood_indicator == 1) {
    mu_est<-rep(0, times=n)
    for(j in 1:s){
      mu_est[site_id == j]<-z[site_id == j,]%*%par0_old[1:2] + 
        x[site_id == j,]%*%Bs_new[j,]+eta_new[j]
    }
    y_hat <- exp(mu_est)
    
    # Generate y_hat
    return(list(
      par0_old=par0_old,
      std_errors=std_errors,
      Y_hat=y_hat,
      Bs_means=Bs_means,
      chain_Bs=chain_Bs,
      Bs_sd=Bs_sd,
      Bs_quantile_025=Bs_quantile_025,
      Bs_quantile_975=Bs_quantile_975
    ))
  }
}



