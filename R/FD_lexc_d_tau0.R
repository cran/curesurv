FD_lexc_d_tau0 <- function(z_tau, z_alpha, x, theta) {
  n_z_tau <- ncol(z_tau)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha <- ncol(z_alpha)
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- theta[1]
  if (n_z_tau > 0) {
    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- (theta[n_z_alpha + 2])
    tau0 <- theta[n_z_tau + 3]
    tau_z <- theta[(n_z_tau + 4):(n_z_tau + 4 + n_z_tau_ad)]
    alpha <- exp(alpha0 + z_alpha %*% alpha_k)
    tau <- exp(tau0 + z_tau %*% tau_z)
    u <- x / (tau)

  }else {
    beta <- (theta[2])
    tau0 <- theta[3]
    alpha <- exp(alpha0)
    tau <- exp(tau0)
    u <- x / (tau)
  }

  res <- ifelse((x <= (tau)),
                -(u ^ ((alpha) - 1) * (1 - u) ^ (exp(beta - 1)))*((alpha-1)+(u/(1-u))*exp(beta - 1)),
                0)
  return(res)
}
