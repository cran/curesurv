#' @import stats
cumLexc_ad2_topred <- function(z_tau, z_alpha, x, theta)
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


    alpha <- (alpha0 + z_alpha %*% alpha_k)
    tau <- (tau0 + z_tau %*% tau_z)
    u <- x/(tau)
    beta2 <- beta
    cumhaz <- ifelse((x <= (tau)),
                     (tau) * beta(alpha, beta2) * pbeta(u, alpha, beta2),
                     (tau) * beta(alpha, beta2) * pbeta(1, alpha, beta2)
    )
    cumhaz2 <- (tau) * beta(alpha, beta2) * pbeta(1, alpha, beta2)

  } else if (n_z_tau == 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- (theta[n_z_alpha + 2])
    tau0 <- theta[n_z_alpha + 2 + 1]
    alpha <- (alpha0 + z_alpha %*% alpha_k)
    tau <- (tau0 )

    u <- x / (tau)

    beta2 <- beta
    cumhaz <- ifelse((x <= (tau)),
                     (tau) * beta(alpha, beta2) * pbeta(u, alpha, beta2),
                     (tau) * beta(alpha, beta2) * pbeta(1, alpha, beta2)
    )
    cumhaz2 <- (tau) * beta(alpha, beta2) * pbeta(1, alpha, beta2)

  }else if (n_z_tau > 0 & n_z_alpha == 0) {

    beta <- (theta[2])
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]

    alpha <- (alpha0)

    tau <- (tau0 + z_tau %*% tau_z)

    u <- x / (tau)
    u <- u*(u < 1)

    beta2 <- beta
    cumhaz <- ifelse((x <= (tau)),
                     (tau) * beta(alpha, beta2) * pbeta(u, alpha, beta2),
                     (tau) * beta(alpha, beta2) * pbeta(1, alpha, beta2)
    )
    cumhaz2 <- (tau) * beta(alpha, beta2) * pbeta(1, alpha, beta2)

  }
  else if (n_z_tau == 0 & n_z_alpha == 0) {
    beta <- (theta[2])
    tau0 <- theta[3]
    alpha <- alpha0
    tau <- (tau0)

    u <- x / (tau)
    u <- u*(u < 1)

    beta2 <- beta
     cumhaz <- ifelse((x <= (tau)),
                      (tau) * beta(alpha, beta2) * pbeta(u, alpha, beta2),
                      (tau) * beta(alpha, beta2) * pbeta(1, alpha, beta2)
     )
     cumhaz2 <- (tau) * beta(alpha, beta2) * pbeta(1, alpha, beta2)
  }
  return(list(netsurv = exp(-cumhaz),
              cumhaz = cumhaz,
              cumhaz2 = cumhaz2,
              pt_cure = (exp(-cumhaz2)/exp(-cumhaz)),
              pi = exp(-cumhaz2),
              tau = tau))
}
