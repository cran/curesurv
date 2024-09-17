#' @title lexc_alphaweibull function
#'
#' @description calculates the instantaneous excess hazard from a Weibull distribution
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#'
#' @param z_pcured covariates matrix acting on cure proportion.
#'
#'
#' @param x the time arguments at which to calculate the cumulative excess hazard
#'
#'
#' @param theta the parameters of the cumulative excess hazard from a
#' Weibull distribution
#'
#' @param  sign_delta only used for mixture cure rate models to specify if the
#' effects or minus the effects of covariates acting on uncured survival to be
#' considered. Default will be sign_delta = "1". The alternative is
#' sign_delta = "-1".
#'
#' @keywords lexc_alphaweibull
#'
#' @return An object of class \code{curesurv}.
#' This object is a vector containing:
#'
#' @references Mudholkar, G.S. and Srivastava, D.K. (1993).
#' Exponentiated Weibull family for analyzing bathtub failure-rate data,
#' IEEE Transactions on Reliability, 42, 299–302.
#'
#'
#'
#' Mudholkar, G.S., Srivastava, D.K., and Freimer, M. (1995). The exponentiated
#' Weibull family: a reanalysis of the bus-motor-failure data,
#'  Technometrics, 37, 436-445.doi:10.2307/1269735
#' (\href{https://www.jstor.org/stable/1269735}{jstor})
#'
#' @keywords internal


lexc_alphaweibull <- function(z_ucured =  z_ucured, z_pcured = z_pcured,
                              x = x, theta = theta, sign_delta = 1) {
  n_z_pcured <- ncol(z_pcured)
  n_z_ucured <- ncol(z_ucured)

  if (n_z_pcured > 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <-sign_delta* theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0 + z_pcured %*% betak
    cured <- 1/(1 + exp(-pcure))
    usurv <- (exp(-exp(lambda)*(x)^exp(gamma)))^exp(z_ucured %*% delta)
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1)) * exp(z_ucured %*% delta)
    u_f <- uhaz*usurv
    SurvE <- cured + (1 - cured)*usurv
    cumHazE <- -log(SurvE)
    ehaz <- ((1 - cured)*u_f) / (cured + (1 - cured)*usurv)

  }else if (n_z_pcured > 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <-sign_delta* theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0 + z_pcured %*% betak
    cured <- 1/(1 + exp(-pcure))
    usurv <- (exp(-exp(lambda)*(x)^exp(gamma)))
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1))
    u_f <- uhaz*usurv
    SurvE <- cured + (1 - cured)*usurv
    cumHazE <- -log(SurvE)
    ehaz <- ((1 - cured)*u_f) / (cured + (1 - cured)*usurv)

  }else if (n_z_pcured == 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <-sign_delta* theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0
    cured <- 1/(1 + exp(-pcure))
    usurv <- (exp(-exp(lambda)*(x)^exp(gamma)))^exp(z_ucured %*% delta)
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1)) * exp(z_ucured %*% delta)
    u_f <- uhaz*usurv
    SurvE <- cured + (1 - cured)*usurv
    cumHazE <- -log(SurvE)
    ehaz <- ((1 - cured)*u_f) / (cured + (1 - cured)*usurv)

  } else if (n_z_pcured == 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    lambda <- theta[2]
    gamma <- theta[3]
    pcure <- beta0
    cured <- 1/(1 + exp(-pcure))
    usurv <- (exp(-exp(lambda)*(x)^exp(gamma)))
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1))
    u_f <- uhaz*usurv
    SurvE <- cured + (1 - cured)*usurv
    cumHazE <- -log(SurvE)
    ehaz <- ((1 - cured)*u_f) / (cured + (1 - cured)*usurv)
  }
  return(ehaz)
}
