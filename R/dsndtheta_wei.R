#' @title dsndtheta_wei function
#'
#' @description function of partial derivates of net survival depending on
#' theta (estimated parameters) from a mixture cure model with uncured survival
#' following a Weibull distribution
#'
#' @param theta estimated parameters from the ouput of a mixture cure model implemented in curesurv
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param x time at which the estimates are predicted
#'
#' @param cumLexctopred pre prediction obtained with cumLexc_alphaweibull_topred
#'
#' @keywords internal

dsndtheta_wei <- function(z_ucured =  z_ucured, z_pcured = z_pcured,
                          x = x, theta = theta,cumLexctopred)
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
    cured <- cumLexctopred$cured
    usurv <- cumLexctopred$usurv
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1)) * exp(z_ucured %*% delta)
    u_f <- uhaz*usurv
    SurvE <- cumLexctopred$SurvE
    cumHazE <- cumLexctopred$cumhaz

    dsndbeta0 <- (cured*(1 - usurv))/(1 + exp(pcure))
    dsndbetak <-  sweep(z_pcured, MARGIN = 1, dsndbeta0, '*')
    dsndlambda <-  -(x^exp(gamma) * exp(z_ucured %*% delta) * usurv * cured  * exp(lambda))/(exp(pcure))
    dsndgamma <- dsndlambda * log(x) * exp(gamma)
    dsnddelta <-  sweep(z_ucured, MARGIN = 1, dsndlambda, '*')

    derivees_partielles <- cbind(dsndbeta0 = dsndbeta0,
                                 dsndbetak = dsndbetak,
                                 dsndlambda = dsndlambda,
                                 dsndgamma = dsndgamma,
                                 dsnddelta = dsnddelta)

  } else if (n_z_pcured > 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]

    pcure <- beta0 + z_pcured %*% betak
    cured <- cumLexctopred$cured
    usurv <- cumLexctopred$usurv
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1))
    u_f <- uhaz*usurv
    SurvE <- cumLexctopred$SurvE
    cumHazE <- cumLexctopred$cumhaz


    dsndbeta0 <- (cured*(1 - usurv))/(1 + exp(pcure))
    dsndbetak <-  sweep(z_pcured, MARGIN = 1, dsndbeta0, '*')
    dsndlambda <-  -(x^exp(gamma) * usurv * cured  * exp(lambda))/(exp(pcure))
    dsndgamma <- dsndlambda * log(x) * exp(gamma)

    derivees_partielles <- cbind(dsndbeta0 = dsndbeta0,
                                 dsndbetak = dsndbetak,
                                 dsndlambda = dsndlambda,
                                 dsndgamma = dsndgamma)


  } else if (n_z_pcured == 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0
    cured <- cumLexctopred$cured
    usurv <- cumLexctopred$usurv
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1)) * exp(z_ucured %*% delta)
    u_f <- uhaz*usurv
    SurvE <- cumLexctopred$SurvE
    cumHazE <- cumLexctopred$cumhaz

    dsndbeta0 <- (cured*(1 - usurv))/(1 + exp(pcure))
    dsndlambda <-  -(x^exp(gamma) * exp(z_ucured %*% delta) * usurv * cured  * exp(lambda))/(exp(pcure))
    dsndgamma <- dsndlambda * log(x) * exp(gamma)
    dsnddelta <-  sweep(z_ucured, MARGIN = 1, dsndlambda, '*')

    derivees_partielles <- cbind(dsndbeta0 = dsndbeta0,
                                 dsndlambda = dsndlambda,
                                 dsndgamma = dsndgamma,
                                 dsnddelta = dsnddelta)






  } else if (n_z_pcured == 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    lambda <- theta[2]
    gamma <- theta[3]
    pcure <- beta0
    cured <- cumLexctopred$cured
    usurv <- cumLexctopred$usurv
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1))
    u_f <- uhaz*usurv
    SurvE <- cumLexctopred$SurvE
    cumHazE <- cumLexctopred$cumhaz




    dsndbeta0 <- (cured * (1 - usurv))/(1 + exp(beta0))
    dsndlambda <- (-(cured * x^(exp(gamma))*usurv)/exp(beta0))*exp(lambda)
    dsndgamma <- exp(lambda) * log(x) * dsndlambda * exp(gamma)

    derivees_partielles <- cbind(dsndbeta0 = dsndbeta0,
                                dsndlambda = dsndlambda,
                                dsndgamma = dsndgamma)
  }




  return(derivees_partielles)
}

