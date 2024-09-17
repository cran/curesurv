

FD_LambdaE_d_tau0 <- function(z_tau, z_alpha, x, theta) {
  n_z_tau <- ncol(z_tau)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha <- ncol(z_alpha)
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- theta[1]
  if(n_z_tau>0){
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

  beta2 <- (exp(beta - 1) + 1)
  uprime <- ifelse((u <= 1 & u >= 0), u, 1)

  resdll_tau0 <- (beta(alpha, beta2)*(tau * stats::pbeta(uprime, alpha, beta2) +
                                        tau * stats::dbeta(tau, shape1 = alpha,
                                                  shape2 = beta2,
                                                  ncp = 0, log = FALSE) * tau))
  return(resdll_tau0)
}
