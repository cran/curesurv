#' @title pt_cure_wei function
#'
#' @description calculates the probability Pi(t) of being
#' cured at a given time t after diagnosis knowing that he/she was alive up to
#' time t. The predictions are based on a mixture cure model with weibull
#' distribution for the survival of uncured patients.
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param x time at which the estimates are predicted
#'
#' @param theta estimated parameters of the cumulative excess hazard from a mixture
#'  model using curesurv and uncured survival following a Weibull distribution
#'
#' @keywords internal


pt_cure_wei <- function(z_pcured = z_pcured, z_ucured = z_ucured, x, theta) {
  n_z_pcured <- ncol(z_pcured)
  n_z_ucured <- ncol(z_ucured)
  if (n_z_pcured > 0 & n_z_ucured > 0) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <-
      theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]

    pcure <- beta0 + z_pcured %*% betak
    cured <- 1 / (1 + exp(-pcure))
    usurv <-
      (exp(-exp(lambda) * (x) ^ exp(gamma))) ^ exp(z_ucured %*% delta)
    uhaz <-
      exp(gamma) * exp(lambda) * ((x) ^ (exp(gamma) - 1)) * exp(z_ucured %*% delta)
    u_f <- uhaz * usurv
    SurvE <- cured + (1 - cured) * usurv
    cumHazE <- -log(SurvE)
  } else if (n_z_pcured > 0 & n_z_ucured == 0)
    {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <-
      -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0 + z_pcured %*% betak
    cured <- 1 / (1 + exp(-pcure))
    usurv <- (exp(-exp(lambda) * (x) ^ exp(gamma)))
    uhaz <- exp(gamma) * exp(lambda) * ((x) ^ (exp(gamma) - 1))
    u_f <- uhaz * usurv
    SurvE <- cured + (1 - cured) * usurv
    cumHazE <- -log(SurvE)

  } else if (n_z_pcured == 0 & n_z_ucured > 0)
    {
    beta0 <- theta[1]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <-
      theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0
    cured <- 1 / (1 + exp(-pcure))
    usurv <-
      (exp(-exp(lambda) * (x) ^ exp(gamma))) ^ exp(z_ucured %*% delta)
    uhaz <-
      exp(gamma) * exp(lambda) * ((x) ^ (exp(gamma) - 1)) * exp(z_ucured %*% delta)
    u_f <- uhaz * usurv
    SurvE <- cured + (1 - cured) * usurv
    cumHazE <- -log(SurvE)
  } else if (n_z_pcured == 0 & n_z_ucured == 0)
    {
    beta0 <- theta[1]
    lambda <- theta[2]
    gamma <- theta[3]
    pcure <- beta0
    cured <- 1 / (1 + exp(-pcure))
    usurv <- (exp(-exp(lambda) * (x) ^ exp(gamma)))
    uhaz <- exp(gamma) * exp(lambda) * ((x) ^ (exp(gamma) - 1))
    u_f <- uhaz * usurv
    SurvE <- cured + (1 - cured) * usurv
    cumHazE <- -log(SurvE)
  }
  ptcure <- cured/SurvE
  return(ptcure)
}
