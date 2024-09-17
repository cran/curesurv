#' @title dTTCdtheta_wei function
#'
#' @description function of partial derivates of time-to-cure (TTC) by theta
#'  (estimated parameters) from a mixture cure model with uncured survival
#' following a Weibull distribution
#'
#' @param theta estimated parameters from the ouput of a mixture cure model implemented in curesurv
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param epsilon value fixed by user to estimate the TTC \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
#'   By default \eqn{\epsilon = 0.05}.
#'
#' @param TTC time-to-cure previously estimated using TTC_wei
#'
#' @keywords internal


dTTCdtheta_wei <- function(z_ucured =  z_ucured,
                           z_pcured = z_pcured,
                           theta = theta,
                           epsilon = epsilon, TTC)
{

  n_z_pcured <- ncol(z_pcured)
  n_z_pcured <- ncol(z_pcured)
  n_z_ucured <- ncol(z_ucured)
  if (n_z_pcured > 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0 + z_pcured %*% betak


    time_to_cure_ttc <- TTC

    dTTCdbeta0 <-
      (1 / exp(gamma)) * time_to_cure_ttc * log((epsilon / (1 - epsilon)) * exp(pcure)) ^
      (-1)
    dTTCdbetak <- sweep(z_pcured, MARGIN = 1, dTTCdbeta0, '*')
    dTTCdlambda <- -(1 / exp(gamma)) * time_to_cure_ttc
    dTTCdgamma <-
      -(1 / exp(gamma)) * log(-(1 / exp(lambda)) * log(((
        epsilon / (1 - epsilon)
      ) * exp(pcure)) ^ exp(-z_ucured %*% delta))) * time_to_cure_ttc
    dTTCddelta <- sweep(z_ucured, MARGIN = 1, dTTCdlambda, '*')

    derivees_partielles <- cbind(
      dTTCdbeta0 = dTTCdbeta0,
      dTTCdbetak = dTTCdbetak,
      dTTCdlambda = dTTCdlambda,
      dTTCdgamma = dTTCdgamma,
      dTTCddelta = dTTCddelta
    )



  } else if (n_z_pcured > 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0 + z_pcured %*% betak


    time_to_cure_ttc <- TTC

    dTTCdbeta0 <-
      (1 / exp(gamma)) * time_to_cure_ttc * log((epsilon / (1 - epsilon)) * exp(pcure)) ^
      (-1)
    dTTCdbetak <- sweep(z_pcured, MARGIN = 1, dTTCdbeta0, '*')
    dTTCdlambda <- -(1 / exp(gamma)) * time_to_cure_ttc
    dTTCdgamma <-
      -(1 / exp(gamma)) * log(-(1 / exp(lambda)) * log(((
        epsilon / (1 - epsilon)
      ) * exp(pcure)) )) * time_to_cure_ttc


    derivees_partielles <- cbind(
      dTTCdbeta0 = dTTCdbeta0,
      dTTCdbetak = dTTCdbetak,
      dTTCdlambda = dTTCdlambda,
      dTTCdgamma = dTTCdgamma
    )


  } else if (n_z_pcured == 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0

    time_to_cure_ttc <- TTC

    dTTCdbeta0 <-
      (1 / exp(gamma)) * time_to_cure_ttc * log((epsilon / (1 - epsilon)) * exp(pcure)) ^(-1)
    dTTCdlambda <- -(1 / exp(gamma)) * time_to_cure_ttc
    dTTCdgamma <- -(1 / exp(gamma)) * log(-(1 / exp(lambda)) * log(((
        epsilon / (1 - epsilon)
      ) * exp(pcure)) ^ exp(-z_ucured %*% delta))) * time_to_cure_ttc

    dTTCddelta <- sweep(z_ucured, MARGIN = 1, dTTCdlambda, '*')

    derivees_partielles <- cbind(
      dTTCdbeta0 = dTTCdbeta0,
      dTTCdlambda = dTTCdlambda,
      dTTCdgamma = dTTCdgamma,
      dTTCddelta = dTTCddelta
    )

  } else if (n_z_pcured == 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    lambda <- theta[2]
    gamma <- theta[3]
    pcure <- beta0
    time_to_cure_ttc <- TTC

    dTTCdbeta0 <-(1 / exp(gamma)) * time_to_cure_ttc * log((epsilon / (1 - epsilon)) * exp(pcure)) ^
      (-1)

    dTTCdlambda <- -(1 / exp(gamma)) * time_to_cure_ttc
    dTTCdgamma <--(1 / exp(gamma)) * log(-(1 / exp(lambda)) * log(((
        epsilon / (1 - epsilon)
      ) * exp(pcure)) ^ exp(0))) * time_to_cure_ttc


    derivees_partielles <- cbind(dTTCdbeta0 = dTTCdbeta0,
                                dTTCdlambda = dTTCdlambda,
                                dTTCdgamma = dTTCdgamma)
  }




  return(derivees_partielles)
}
