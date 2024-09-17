
lexc_mul <- function(z_tau, z_alpha, x, theta) {
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
    alpha <- (exp((alpha0) + z_alpha %*% alpha_k))
    tau <- (exp(tau0 + z_tau %*% tau_z))
    u <- x / (tau)
    haz <- ifelse((x >= 0 & x <= (tau)),
                  (u ^ ((alpha) - 1) * (1 - u) ^ ((beta - 1))),
                  0)
  } else if (n_z_tau > 0 & n_z_alpha == 0) {
    beta <- exp(theta[2]) + 1
    tau0 <- theta[2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    alpha <- exp(alpha0)
    tau <- (exp(tau0 + z_tau %*% tau_z))
    u <- x / (tau)
    haz <- ifelse((x >= 0 & x <= (tau)),
                  (u ^ ((alpha) - 1) * (1 - u) ^ ((beta - 1))),
                  0)
  } else if (n_z_tau == 0 & n_z_alpha > 0) {
    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- exp(theta[n_z_alpha + 2]) + 1
    tau0 <- theta[n_z_alpha + 2 + 1]
    alpha <- (exp(alpha0 + z_alpha %*% alpha_k))
    tau <- tau0
    u <- x / (tau)
    haz <- ifelse((x >= 0 & x <= (tau)),
                  (u ^ ((alpha) - 1) * (1 - u) ^ ((beta - 1))),
                  0)
  } else if (n_z_tau == 0 & n_z_alpha == 0)  {
    beta <- exp(theta[2]) + 1
    tau0 <- theta[3]
    alpha <- exp(alpha0)
    tau <- exp(tau0)
    u <- x / (tau)
    haz <- ifelse((x <= (tau)),
                  (u ^ ((alpha) - 1) * (1 - u) ^ ((beta - 1))),
                  0)

  }
  return(haz)
}


