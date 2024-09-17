#' @title cumLexc_mul function
#'
#' @description returns the cumulative excess hazard for an TNEH model
#'  in case of parametrization of log the of the time to null excess hazard as
#'  function to fit the data
#'
#' @param z_tau covariates depending on tau
#' @param z_alpha covariates depending on alpha
#' @param x time value
#'
#' @param theta of the coefficient of tneh parameters
#'
#' @keywords cumLexc_mul
#'
#' @return An object of class numeric containing the cumulative excess
#' hazard with the same length as the time.
#'
#'
#' @export
#'
cumLexc_mul <- function(z_tau, z_alpha, x, theta)
{
  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- theta[1]
  if (n_z_tau > 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- theta[n_z_alpha + 2]
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    alpha <- (exp(alpha0 + z_alpha %*% alpha_k))
    tau <- (exp(tau0 + z_tau %*% tau_z))
    u <- x / (tau)
    beta2 <- exp(beta) + 1
    cumhaz <- ifelse((x <= (tau)),
                     (tau) * beta(alpha,
                                  beta2) *
                       stats::pbeta(u, alpha, beta2),
                     (tau) * beta(alpha,
                                  beta2) *
                       stats::pbeta(1, alpha, beta2)
    )
  } else if (n_z_tau > 0 & n_z_alpha == 0) {

    beta <- theta[n_z_alpha + 2]
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    alpha <- exp(alpha0)
    tau <- (exp(tau0 + z_tau %*% tau_z))
    u <- x / (tau)
    beta2 <- exp(beta) + 1
    cumhaz <- ifelse((x <= (tau)),
                     (tau) * beta(alpha,
                                  beta2) *
                       stats::pbeta(u, alpha, beta2),
                     (tau) * beta(alpha,
                                  beta2) *
                       stats::pbeta(1, alpha, beta2)
    )

  } else if (n_z_tau == 0 & n_z_alpha > 0) {
    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- theta[n_z_alpha + 2]
    tau0 <- theta[n_z_alpha + 2 + 1]
    alpha <- (exp(alpha0 + z_alpha %*% alpha_k))
    tau <- (exp(tau0))
    u <- x / (tau)
    beta2 <- exp(beta) + 1
    cumhaz <- ifelse((x <= (tau)),
                     (tau) * beta(alpha,
                                  beta2) *
                       stats::pbeta(u, alpha, beta2),
                     (tau) * beta(alpha,
                                  beta2) *
                       stats::pbeta(1, alpha, beta2)
    )

  } else if (n_z_tau == 0 & n_z_alpha == 0) {
    beta <- (theta[2])
    tau0 <- theta[3]
    alpha <- exp(alpha0)
    tau <- exp(tau0)
    u <- x / (tau)
    beta2 <- (exp(beta) + 1)
    cumhaz <- ifelse((x <= (tau)),
                     (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                     (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
    )
  }
  return(cumhaz)
}

