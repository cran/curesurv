cumLexc_m1_mul <- function(z_tau, z_alpha,z_c, x, theta)
{
  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1
  n_z_c <- ncol(z_c)
  n_z_c_ad <- n_z_c - 1
  alpha0 <- theta[1]
  if (n_z_tau > 0 & n_z_alpha > 0) {
    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- (theta[n_z_alpha + 2])
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    if (n_z_c > 0) {
      c_par0 <- theta[(n_z_alpha + 2 + n_z_tau + 1 + 1)]
      c_k <- theta[(n_z_alpha + 2 + n_z_tau + 1 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1 + 1 + n_z_c)]
      c_par <- rowSums(c_par0 + z_c %*% c_k)
    }else{
      c_par <- c_par0 <- theta[(n_z_alpha + 2 + n_z_tau + 1 + 1)]
    }
    alpha <- rowSums(exp(alpha0 + z_alpha %*% alpha_k))
    tau <- rowSums(exp(tau0 + z_tau %*% tau_z))
    u <- x / (tau)
    beta2 <- (exp(beta - 1) + 1)
    cumhaz <- ifelse((x <= (tau)),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2) + exp(c_par) * x),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2) + exp(c_par) * x)
    )
  } else if (n_z_tau > 0 & n_z_alpha == 0) {
    beta <- (theta[n_z_alpha + 2])
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    if (n_z_c > 0) {
      c_par0 <- theta[(n_z_alpha + 2 + n_z_tau + 1 + 1)]
      c_k <- theta[(n_z_alpha + 2 + n_z_tau + 1 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1 + 1 + n_z_c)]
      c_par <- rowSums(c_par0 + z_c %*% c_k)
    }else{
      c_par <- c_par0 <- theta[(n_z_alpha + 2 + n_z_tau + 1 + 1)]
    }
    alpha <- exp(alpha0)
    tau <- rowSums(exp(tau0 + z_tau %*% tau_z))
    u <- x / (tau)
    beta2 <- (exp(beta - 1) + 1)
    cumhaz <- ifelse((x <= (tau)),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2) + exp(c_par) * x),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2) + exp(c_par) * x)
    )
  } else if (n_z_tau == 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- (theta[n_z_alpha + 2])
    tau0 <- theta[n_z_alpha + 2 + 1]
    if (n_z_c > 0) {
      c_par0 <- theta[(n_z_alpha + 2 + n_z_tau + 1 + 1)]
      c_k <- theta[(n_z_alpha + 2 + n_z_tau + 1 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1 + 1 + n_z_c)]
      c_par <- rowSums(c_par0 + z_c %*% c_k)
    }else{
      c_par <- c_par0 <- theta[(n_z_alpha + 2 + n_z_tau + 1 + 1)]
    }
    alpha <- rowSums(exp(alpha0 + z_alpha %*% alpha_k))
    tau <- exp(tau0)
    u <- x / (tau)
    beta2 <- (exp(beta - 1) + 1)
    cumhaz <- ifelse((x <= (tau)),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2) + exp(c_par) * x),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2) + exp(c_par) * x)
    )

  } else {
    beta <- (theta[2])
    tau0 <- theta[3]
    c_par <- theta[4]
    alpha <- exp(alpha0)
    tau <- exp(tau0)
    u <- x / (tau)
    beta2 <- (exp(beta - 1) + 1)

    cumhaz <- ifelse((x <= (tau)),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2) + exp(c_par) * x),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2) + (exp(c_par) * x))
    )
  }
  return(cumhaz)
}
