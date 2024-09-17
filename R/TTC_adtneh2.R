#' @title TTC_adtneh2 function
#'
#' @description calculates the probability Pi(t) of being cured at a given time
#' t after diagnosis knowing that he/she was alive up to time t from a
#' Time-to-Null excess hazard model using numerical method, uniroot. In other
#' words,
#' Pi(t)=(probability of being cured and alive up to time t given xi)/
#'  (probability of being alive up to time t given xi)
#'
#' @param object ouput from a non mixture model with distribution "tneh"
#' from curesurv function.
#'
#' @param z_alpha Covariates matrix acting on parameter alpha of the density of
#'  time-to-null excess hazard model
#'
#' @param z_tau Covariates matrix acting on time-to-null parameter.
#'
#' @param xmax time max at which Pi(t) is calculated.
#'
#' @param epsilon value fixed by user to estimate the TTC \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
#'   By default \eqn{\epsilon = 0.05}.
#'
#' @author Juste Goungounga, Judith Breaud Eugenie Blandin, Olayide Boussari, Valerie Jooste
#'
#' @references Boussari O, Bordes L, Romain G, Colonna M, Bossard N, Remontet L,
#'  Jooste V. Modeling excess hazard with time-to-cure as a parameter.
#'   Biometrics. 2021 Dec;77(4):1289-1302. doi: 10.1111/biom.13361.
#'    Epub 2020 Sep 12. PMID: 32869288.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/32869288/}{pubmed})
#'
#'
#'  Boussari O, Romain G, Remontet L, Bossard N, Mounier M, Bouvier AM,
#'  Binquet C, Colonna M, Jooste V. A new approach to estimate time-to-cure from
#'  cancer registries data. Cancer Epidemiol. 2018 Apr;53:72-80.
#'  doi: 10.1016/j.canep.2018.01.013. Epub 2018 Feb 4. PMID: 29414635.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/29414635/}{pubmed})
#'
#' @keywords internal



TTC_adtneh2 <- function(z_alpha, z_tau, xmax, object, epsilon = epsilon) {

  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")

  pt_cure_tneh <- function(z_alpha, z_tau, x, object, epsilon = epsilon) {
    theta <- object$coefficients
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
      u <- x / (tau)
      beta2 <- beta
      cumhaz <- ifelse((x <= (tau)),
                       (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                       (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
      )
      cumhaz2 <- (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)

    } else if (n_z_tau == 0 & n_z_alpha > 0) {

      alpha_k <- theta[2:(n_z_alpha + 1)]
      beta <- (theta[n_z_alpha + 2])
      tau0 <- theta[n_z_alpha + 2 + 1]
      alpha <- (alpha0 + z_alpha %*% alpha_k)
      tau <- (tau0 )

      u <- x / (tau)

      beta2 <- beta
      cumhaz <- ifelse((x <= (tau)),
                       (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                       (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
      )
      cumhaz2 <- (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)

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
                       (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                       (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
      )
      cumhaz2 <- (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)

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
                       (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                       (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
      )
      cumhaz2 <- (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
    }


    fx <- (exp(-cumhaz2) / exp(-cumhaz)) - (1 - epsilon)
    return(fx)
  }



  res_ttc_tneh2 <- do.call("rbind",
                           ttc_tneh2(object = object,
                                     fx = pt_cure_tneh,
                                     xmax = xmax,
                                     z_tau = z_tau,
                                     z_alpha = z_alpha,
                                     epsilon = epsilon))


res_all <- list(TTC = unlist(res_ttc_tneh2[, "root"]),
                f.root = unlist(res_ttc_tneh2[, "f.root"]),
                iter = unlist(res_ttc_tneh2[, "iter"]),
                init.it = unlist(res_ttc_tneh2[, "init.it"]))

  return(res_all)


}
