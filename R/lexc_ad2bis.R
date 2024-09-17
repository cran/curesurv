lexc_ad2bis <- function(z_tau, z_alpha, x, theta) {
  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- theta[1]
  if (n_z_tau > 0 & n_z_alpha > 0) {
    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- exp(theta[n_z_alpha + 2]) + 1
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    alpha <- exp(alpha0 + z_alpha %*% alpha_k)
    tau <- exp(tau0 + z_tau %*% tau_z)
  } else if (n_z_tau > 0 & n_z_alpha == 0)
  {
    beta <- exp(theta[2]) + 1
    tau0 <- theta[2 + 1]
    tau_z <- theta[(2 + 1 + 1):(2 + n_z_tau + 1)]
    alpha <- exp(alpha0) * rep(1,length(x))
    tau <- exp(tau0 + z_tau %*% tau_z)
  } else if (n_z_tau == 0 & n_z_alpha > 0) {
    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- exp(theta[n_z_alpha + 2]) + 1
    tau0 <- theta[n_z_alpha + 2 + 1]
    alpha <- exp(alpha0 + z_alpha %*% alpha_k)
    tau <- exp(tau0) * rep(1,length(x))
  } else
    if (n_z_tau == 0 & n_z_alpha == 0) {

      beta <- exp(theta[2]) + 1
      tau0 <- theta[3]
      alpha <- exp(alpha0) * rep(1,length(x))
      tau <- exp(tau0) * rep(1,length(x))
    }

  u <- x / (tau)
  haz <- rep(0, length(x))
  haz[which(x <= tau)] <- (u[which(x <= tau)]) ^ ((alpha[which(x <= tau)]) - 1) * (1 - u[which(x <= tau)]) ^ (beta - 1)
  return(haz)
}

