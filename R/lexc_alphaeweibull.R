#' @title lexc_alphaeweibull function
#'
#' @description calculates the instantaneous excess hazard from an exponentiated Weibull distribution
#'
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
#' @param theta the parameters of the cumulative excess hazard from an exponentiated Weibull distribution
#'
#' @param  sign_delta only used for mixture cure rate models to specify if the
#' effects or minus the effects of covariates acting on uncured survival to be
#' considered. Default will be sign_delta = "1". The alternative is
#' sign_delta = "-1".
#'
#'
#' @keywords lexc_alphaeweibull
#'
#' @return An object of class \code{curesurv}.
#' This object is a vector containing:
#'
#' @references Mudholkar, G.S. and Srivastava, D.K. (1993).
#' Exponentiated Weibull family for analyzing bathtub failure-rate data,
#' IEEE Transactions on Reliability, 42, 299-302.
#'
#'
#'
#' Mudholkar, G.S., Srivastava, D.K., and Freimer, M. (1995). The exponentiated
#' Weibull family: a reanalysis of the bus-motor-failure data,
#'  Technometrics, 37, 436–445.doi:10.2307/1269735
#' (\href{https://www.jstor.org/stable/1269735}{jstor})
#'
#'
#' @keywords internal


lexc_alphaeweibull <- function(z_ucured =  z_ucured, z_pcured = z_pcured,
                               x = x, theta = theta,  sign_delta = 1) {
  n_z_pcured <- ncol(z_pcured)
  n_z_ucured <- ncol(z_ucured)
  if (n_z_pcured > 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    thetapw <- theta[(1 + n_z_pcured + 3)]
    delta <-sign_delta*theta[(1 + n_z_pcured + 4):(1 + n_z_pcured + 3 + n_z_ucured)]
    pcure <- beta0 + z_pcured %*% betak
    cured <- 1/(1 + exp(-pcure))
    usurv <- (1 - (1 - exp(-(exp(lambda+z_ucured %*% delta) * x) ^ exp(gamma))) ^ exp(thetapw))
    u_f <-exp(gamma) * exp(thetapw) * exp(lambda+z_ucured %*% delta) *
      (1 - exp(-(exp(lambda+z_ucured %*% delta) * x ) ^ exp(gamma))) ^ (exp(thetapw) - 1) *
      exp(-(exp(lambda+z_ucured %*% delta) * x) ^ exp(gamma)) *
      (exp(lambda+z_ucured %*% delta)*x) ^ (exp(gamma) - 1)
    uhaz <- u_f/usurv
    SurvE <- cured + (1 - cured)*usurv
    cumHazE <- -log(SurvE)
    ehaz <- ((1 - cured)*u_f) / (cured + (1 - cured)*usurv)

  } else if (n_z_pcured == 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    thetapw <- theta[(1 + n_z_pcured + 3)]
    delta <-sign_delta*theta[(1 + n_z_pcured + 4):(1 + n_z_pcured + 3 + n_z_ucured)]
    pcure <- beta0
    cured <- 1/(1 + exp(-pcure))
    usurv <- (1 - (1 - exp(-(exp(lambda+z_ucured %*% delta) * x) ^ exp(gamma))) ^ exp(thetapw))
    u_f <-exp(gamma) * exp(thetapw) * exp(lambda+z_ucured %*% delta) *
      (1 - exp(-(exp(lambda+z_ucured %*% delta) * x ) ^ exp(gamma))) ^ (exp(thetapw) - 1) *
      exp(-(exp(lambda+z_ucured %*% delta) * x) ^ exp(gamma)) *
      (exp(lambda+z_ucured %*% delta)*x) ^ (exp(gamma) - 1)
    uhaz <- u_f/usurv
    SurvE <- cured + (1 - cured)*usurv
    cumHazE <- -log(SurvE)
    ehaz <- ((1 - cured)*u_f) / (cured + (1 - cured)*usurv)

  } else if (n_z_pcured > 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    thetapw <- theta[(1 + n_z_pcured + 3)]
    pcure <- beta0 + z_pcured %*% betak
    cured <- 1/(1 + exp(-pcure))
    usurv <- (1 - (1 - exp(-(exp(lambda) * x) ^ exp(gamma))) ^ exp(thetapw))
    u_f <-exp(gamma) * exp(thetapw) * exp(lambda) *
      (1 - exp(-(exp(lambda) * x ) ^ exp(gamma))) ^ (exp(thetapw) - 1) *
      exp(-(exp(lambda) * x) ^ exp(gamma)) *
      (exp(lambda)*x) ^ (exp(gamma) - 1)

    uhaz <- u_f/usurv
    SurvE <- cured + (1 - cured)*usurv
    cumHazE <- -log(SurvE)
    ehaz <- ((1 - cured)*u_f) / (cured + (1 - cured)*usurv)
  } else if (n_z_pcured == 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    lambda <- theta[2]
    gamma <- theta[3]
    thetapw <- theta[4]
    pcure <- beta0
    cured <- 1/(1 + exp(-pcure))
    usurv <- (1 - (1 - exp(-(exp(lambda) * x)^exp(gamma)))^exp(thetapw))

    u_f <-exp(gamma) * exp(thetapw) * exp(lambda) *
      (1 - exp(-(exp(lambda) * x ) ^ exp(gamma))) ^ (exp(thetapw) - 1) *
      exp(-(exp(lambda) * x) ^ exp(gamma)) *
      (exp(lambda)*x) ^ (exp(gamma) - 1)


    uhaz <- u_f/usurv
    SurvE <- cured + (1 - cured)*usurv
    cumHazE <- -log(SurvE)
    ehaz <- ((1 - cured)*u_f) / (cured + (1 - cured)*usurv)
  }
  return(ehaz)
}
