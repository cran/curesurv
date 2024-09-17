
cumLexc_m1_ad_topred <- function(z_tau, z_alpha,z_c, x, theta)
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
      c_par0 <- (theta[(n_z_alpha + 2 + n_z_tau + 1 + 1)])
      c_k <- theta[(n_z_alpha + 2 + n_z_tau + 1 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1 + 1 + n_z_c)]
      c_par <- exp(rowSums(c_par0 + z_c %*% c_k))
    }else{
      c_par <- c_par0 <- exp(theta[(n_z_alpha + 2 + n_z_tau + 1 + 1)])
    }
    alpha <- abs(rowSums(alpha0 + z_alpha %*% alpha_k))
    tau <- abs(rowSums(tau0 + z_tau %*% tau_z))
    u <- x / (tau)
    beta2 <- ((beta - 1) + 1)
    cumhaz <- ifelse((x <= (tau)),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2) + (c_par) * x),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2) + (c_par) * x)
    )
    cumhaz2 <-       ((tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2) + (c_par) * tau)

  } else if (n_z_tau > 0 & n_z_alpha == 0) {
    beta <- (theta[n_z_alpha + 2])
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    if (n_z_c > 0) {
      c_par0 <- (theta[(n_z_alpha + 2 + n_z_tau + 1 + 1)])
      c_k <- theta[(n_z_alpha + 2 + n_z_tau + 1 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1 + 1 + n_z_c)]
      c_par <- exp(rowSums(c_par0 + z_c %*% c_k))
    }else{
      c_par <- c_par0 <- exp(theta[(n_z_alpha + 2 + n_z_tau + 1 + 1)])
    }
    alpha <- abs(alpha0)
    tau <- abs(rowSums(tau0 + z_tau %*% tau_z))
    u <- x / (tau)
    beta2 <- beta
    cumhaz <- ifelse((x <= (tau)),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2) + (c_par) * x),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2) + (c_par) * x)
    )
    cumhaz2 <-       ((tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2) + (c_par) * tau)

  } else if (n_z_tau == 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- (theta[n_z_alpha + 2])
    tau0 <- theta[n_z_alpha + 2 + 1]
    if (n_z_c > 0) {
      c_par0 <- (theta[(n_z_alpha + 2 + n_z_tau + 1 + 1)])
      c_k <- theta[(n_z_alpha + 2 + n_z_tau + 1 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1 + 1 + n_z_c)]
      c_par <- exp(rowSums(c_par0 + z_c %*% c_k))
    }else{
      c_par <- c_par0 <- exp(theta[(n_z_alpha + 2 + n_z_tau + 1 + 1)])
    }
    alpha <- abs(rowSums(alpha0 + z_alpha %*% alpha_k))
    tau <- abs(tau0)
    u <- x / (tau)
    beta2 <- beta
    cumhaz <- ifelse((x <= (tau)),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2) + (c_par) * x),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2) + (c_par) * x)
    )
    cumhaz2 <-       ((tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2) + (c_par) * tau)


  } else {
    beta <- (theta[2])
    tau0 <- theta[3]
    c_par <- theta[4]
    alpha <- abs(alpha0)
    tau <- (tau0)
    u <- x / (tau)
    beta2 <- ((beta - 1) + 1)

    cumhaz <- ifelse((x <= (tau)),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2) + (c_par) * x),
                     ((tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2) + ((c_par) * x))
    )
    cumhaz2 <-       ((tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2) + (c_par) * tau)

  }
  return(list(cumhaz = cumhaz,
              cumhaz2 = cumhaz2,
              tau = tau))
}






