# Calculates covariance matrix using squared-exponential kernel
get_cov_SE <- function(x, ell) {
  # ell is either a scalar or a vector with as many elements as x has columns
  if(length(ell) == 1) {
    ell <- rep(ell, ncol(x))
  } else {
    stopifnot(length(ell) == ncol(x))
  }
  ell_mat <- matrix(ell, nrow(x), length(ell), byrow = TRUE)
  D <- as.matrix(dist(x / ell_mat))
  K <- exp(-0.5 * D^2)
  return(K)
}

# Calculates log marginal likelihood of Gaussian Process
#   Defaults to noise-free observations of y (sigma_sq = 0)
get_log_ml <- function(Xtrain, y, length_scale, amplitude, sigma_sq = 0) {
  stopifnot(amplitude > 0 & length_scale > 0)
  n1 <- nrow(Xtrain)
  K11 <- amplitude * get_cov_SE(Xtrain, ell = length_scale)
  R <- chol(K11 + diag(sigma_sq, n1))
  b <- backsolve(R, forwardsolve(t(R), y))
  log_ml <- -0.5 * sum(y * b) - sum(log(diag(R))) - (n1 / 2) * log(2 * pi)
  return(log_ml)
}

# Optimize log marginal likelihood
#   Optimizes over amplitude and length scale parameters at a fixed sigma_sq (default zero)
get_log_ml_opt <- function(Xtrain, y, sigma_sq = 0) {
  f <- function(x) {
    theta <- exp(x[1])
    ell <- exp(x[-1])
    logml <- tryCatch(get_log_ml(Xtrain, y, amplitude = theta, length_scale = ell),
                      error = function(e) -Inf, warning = function(e) -Inf)
    return(-logml)
  }
  opt <- optim(c(0, rep(0, ncol(Xtrain))), f)
  out <- opt
  out$theta <- exp(opt$par[1])
  out$ell <- exp(opt$par[-1])
  out$par <- NULL
  return(out)
}


# Calculate posterior mean and variance of Gaussian Process at Xtest.
#   Defaults to noise-free observations of y (sigma_sq = 0)
get_GP_predictions <- function(Xtrain, Xtest, y, length_scale = 1, amplitude = 1, sigma_sq = 0) {
  stopifnot(amplitude > 0 & length_scale > 0)
  K <- amplitude * get_cov_SE(rbind(Xtrain, Xtest), ell = length_scale)
  n1 <- nrow(Xtrain)
  n2 <- nrow(Xtest)
  n <- n1 + n2
  K11 <- K[1:n1, 1:n1]
  K22 <- K[(n1 + 1):n, (n1 + 1):n]
  K12 <- K[1:n1, (n1 + 1):n]

  R <- chol(K11 + diag(sigma_sq, n1))
  b <- backsolve(R, forwardsolve(t(R), y))
  post_mean <- t(K12) %*% b
  C <- forwardsolve(t(R), K12)
  post_var <- K22 - t(C) %*% C

  # n1 is the training sample size
  log_ml <- -0.5 * sum(y * b) - sum(log(diag(R))) - (n1 / 2) * log(2 * pi)

  out <- list(post_mean = post_mean,
              post_var = post_var,
              log_ml = log_ml,
              Xtrain = Xtrain,
              Xtest = Xtest,
              y = y,
              length_scale = length_scale,
              sigma_sq = sigma_sq,
              b = b,
              R = R)

  return(out)
}
