
cumLexc_ad_topred <- function(z_tau, z_alpha, x, theta)
{
  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- abs(theta[1])
  if (n_z_tau > 0) {
    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- (theta[n_z_alpha + 2])
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    alpha <- abs(alpha0 + z_alpha %*% alpha_k)
    tau <- (tau0 + z_tau %*% tau_z)
    u <- x / (tau)
    beta2 <- beta
    cumhaz <- ifelse((x <= (tau)),
                     (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                     (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
    )
    cumhaz2 <-       (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)

  } else{
    beta <- (theta[2])
    tau0 <- theta[3]
    alpha <- abs(alpha0)
    tau <- (tau0)
    u <- x / (tau)
    beta2 <- ((beta - 1) + 1)
    cumhaz <- ifelse((x <= (tau)),
                     (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                     (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
    )
    cumhaz2 <-       (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)

  }
  return(list(cumhaz = cumhaz, cumhaz2 = cumhaz2, tau = tau))
}
