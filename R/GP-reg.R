# Calculates covariance matrix using squared-exponential kernel
get_cov_SE <- function(x, ell) {
  D <- as.matrix(dist(x))
  K <- exp(-0.5 * D^2 / ell^2)
  return(K)
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

  log_ml <- -0.5 * sum(y * b) - sum(log(diag(R))) - (n / 2) * log(2 * pi)

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
