
cumLexc_ad3bis <- function(z_tau, z_alpha, x, theta)
{
  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- (theta[1])
  if (n_z_tau > 0 & n_z_alpha > 0) {
    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- (theta[n_z_alpha + 2])
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    alpha <- exp(alpha0 + z_alpha %*% alpha_k)
    tau <- exp(tau0 + z_tau %*% tau_z)
  } else if (n_z_tau == 0 & n_z_alpha > 0) {
    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- (theta[n_z_alpha + 2])
    tau0 <- theta[n_z_alpha + 2 + 1]
    alpha <- exp(alpha0 + z_alpha %*% alpha_k)
    tau <- exp(tau0) * rep(1,length(x))
  }else if (n_z_tau > 0 & n_z_alpha == 0) {
    beta <- (theta[2])
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    alpha <- exp(alpha0)
    tau <- exp(tau0 + z_tau %*% tau_z)
  }
  else if (n_z_tau == 0 & n_z_alpha == 0) {
    beta <- (theta[2])
    tau0 <- theta[3]
    alpha <- exp(alpha0) * rep(1,length(x))
    tau <- exp(tau0) * rep(1,length(x))
    lambda <- theta[4]
    gamma <- theta[5]
    p <- theta[6]
  }

  beta2 <- exp(beta) + 1
  u <- x / (tau)
  cumhaz <- rep(0, length(x))
  cumhaz[which(x <= tau)] <- (1 - p) * ((tau)[which(x <= tau)] * beta(alpha[which(x <= tau)], beta2) * stats::pbeta(u[which(x <= tau)], alpha[which(x <= tau)], beta2)) -
    (log( exp(-exp(lambda)*((x[which(x <= tau)])^exp(gamma))))) * p

  cumhaz[which(x > tau)] <- (1 - p) * ((tau)[which(x > tau)] * beta(alpha[which(x > tau)], beta2) * stats::pbeta(1, alpha[which(x <= tau)], beta2) ) -
    (log( exp(-exp(lambda)*((x[which(x > tau)])^exp(gamma)))))*p


  return(cumhaz)
}
