#' @title dlogTTCdtheta_wei function
#'
#' @description function of partial derivates of log of time-to-cure \code{log(TTC)}
#' by theta (estimated parameters) from a mixture cure model with uncured survival
#' following a Weibull distribution
#'
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param object ouput from a model implemented in curesurv
#'
#' @param epsilon value fixed by user to estimate the TTC \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
#'   By default \eqn{\epsilon = 0.05}.
#'
#' @param TTC time to cure previsouly estimated by TTC_wei
#'
#' @keywords internal

dlogTTCdtheta_wei <- function(z_ucured =  z_ucured,
                              z_pcured = z_pcured,
                              object = object,
                              epsilon = epsilon,TTC) {
  n_z_pcured <- ncol(z_pcured)
  n_z_pcured <- ncol(z_pcured)
  n_z_ucured <- ncol(z_ucured)
  theta <- object$coefficients
  if (n_z_pcured > 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    id_delta <- (1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)

    delta <- -theta[id_delta]
    pcure <- beta0 + z_pcured %*% betak


    time_to_cure_ttc <- TTC


dTTCdtheta <- dTTCdtheta_wei(z_ucured =  z_ucured,
                             z_pcured = z_pcured,
                             theta = object$coefficients,
                             epsilon = epsilon,TTC)

dlogTTDdbeta0 <- (1/time_to_cure_ttc) * dTTCdtheta[,1]
dlogTTCdbetak <- sweep(as.matrix(dTTCdtheta[,2:(1 + n_z_pcured)]), MARGIN = 1, 1/time_to_cure_ttc, '*')
dlogTTCdlambda <- (1/time_to_cure_ttc) * dTTCdtheta[,(1 + n_z_pcured + 1)]
dlogTTCdgamma <- (1/time_to_cure_ttc) * dTTCdtheta[,(1 + n_z_pcured + 2)]
dlogTTCddelta <-
  sweep(as.matrix(dTTCdtheta[, id_delta]), MARGIN = 1, 1 /
          time_to_cure_ttc, '*')

derivees_partielles <- cbind(dlogTTDdbeta0 = dlogTTDdbeta0,
                             dlogTTCdbetak = dlogTTCdbetak,
                             dlogTTCdlambda = dlogTTCdlambda,
                             dlogTTCdgamma = dlogTTCdgamma,
                             dlogTTCddelta = dlogTTCddelta)



  } else if (n_z_pcured > 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    pcure <- beta0 + z_pcured %*% betak


    time_to_cure_ttc <- TTC


dTTCdtheta <- dTTCdtheta_wei(z_ucured =  z_ucured,
                             z_pcured = z_pcured,
                             theta = object$coefficients,
                             epsilon = epsilon,TTC)

dlogTTDdbeta0 <- (1/time_to_cure_ttc) * dTTCdtheta[,1]
dlogTTCdbetak <- sweep(as.matrix(dTTCdtheta[,2:(1 + n_z_pcured)]),
                       MARGIN = 1, 1/time_to_cure_ttc, '*')
dlogTTCdlambda <- (1/time_to_cure_ttc) * dTTCdtheta[,(1 + n_z_pcured + 1)]
dlogTTCdgamma <- (1/time_to_cure_ttc) * dTTCdtheta[,(1 + n_z_pcured + 2)]


derivees_partielles <- cbind(dlogTTDdbeta0 = dlogTTDdbeta0,
                             dlogTTCdbetak = dlogTTCdbetak,
                             dlogTTCdlambda = dlogTTCdlambda,
                             dlogTTCdgamma = dlogTTCdgamma)


  } else if (n_z_pcured == 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    id_delta <- (1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)
    delta <- -theta[id_delta]
    pcure <- beta0

    time_to_cure_ttc <- TTC


    dTTCdtheta <- dTTCdtheta_wei(z_ucured =  z_ucured,
                                 z_pcured = z_pcured,
                                 theta = object$coefficients,
                                 epsilon = epsilon,TTC)

    dlogTTDdbeta0 <- (1/time_to_cure_ttc) * dTTCdtheta[,1]
    dlogTTCdlambda <- (1/time_to_cure_ttc) * dTTCdtheta[,(1 + n_z_pcured + 1)]
    dlogTTCdgamma <- (1/time_to_cure_ttc) * dTTCdtheta[,(1 + n_z_pcured + 2)]
    dlogTTCddelta <-
      sweep(as.matrix(dTTCdtheta[, id_delta]), MARGIN = 1, 1 /
              time_to_cure_ttc, '*')

    derivees_partielles <- cbind(dlogTTDdbeta0 = dlogTTDdbeta0,
                                 dlogTTCdlambda = dlogTTCdlambda,
                                 dlogTTCdgamma = dlogTTCdgamma,
                                 dlogTTCddelta = dlogTTCddelta)


  } else if (n_z_pcured == 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    lambda <- theta[2]
    gamma <- theta[3]
    pcure <- beta0
    time_to_cure_ttc <- TTC

    dTTCdtheta <- dTTCdtheta_wei(z_ucured =  z_ucured,
                                 z_pcured = z_pcured,
                                 theta = object$coefficients,
                                 epsilon = epsilon,TTC)

    dlogTTCdbeta0 <- (1/time_to_cure_ttc) * dTTCdtheta[,1]
    dlogTTCdlambda <- (1/time_to_cure_ttc) * dTTCdtheta[,2]
    dlogTTCdgamma <- (1/time_to_cure_ttc) * dTTCdtheta[,3]


    derivees_partielles <- cbind(dlogTTCdbeta0 = dlogTTCdbeta0,
                                dlogTTCdlambda = dlogTTCdlambda,
                                dlogTTCdgamma = dlogTTCdgamma)
  }


  return(derivees_partielles)
  }


