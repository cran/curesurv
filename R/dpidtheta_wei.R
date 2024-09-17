#' @title dpitheta_wei function
#'
#' @description Produce partial derivatives of pi the cure proportion
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param x time at which the estimates are predicted
#'
#' @param theta estimated parameters from a mixture model using curesurv
#' and uncured survival following a Weibull distribution
#'
#' @param cumLexctopred description
#'
#' @keywords internal

dpidtheta_wei <- function(z_pcured = z_pcured,
                         z_ucured = z_ucured,
                         x = x, theta,cumLexctopred) {

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

    dpdbeta0 <- -exp(-pcure)*(cured)^2

    dpdbeta_star <- matrix(rep(dpdbeta0,n_z_pcured),ncol=n_z_pcured) * z_pcured


    derivees_partielles <- list(dpdbeta0 = dpdbeta0,
                                dpdbeta_star = dpdbeta_star)

  } else if (n_z_pcured > 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0 + z_pcured %*% betak
    cured <- cumLexctopred$cured
    usurv <- cumLexctopred$usurv
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1))
    u_f <- uhaz*usurv
    SurvE <- cumLexctopred$SurvE
    cumHazE <- cumLexctopred$cumhaz


    dpdbeta0 <- -exp(-pcure)*(cured)^2

    dpdbeta_star <- matrix(rep(dpdbeta0,n_z_pcured),ncol=n_z_pcured) * z_pcured

    derivees_partielles <- list(dpdbeta0 = dpdbeta0,
                                dpdbeta_star = dpdbeta_star)

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

    dpdbeta0 <- -exp(-pcure)*(cured)^2



    derivees_partielles <- list(dpdbeta0 = dpdbeta0)




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


    dpdbeta0 <- -exp(-pcure)*(cured)^2

    derivees_partielles <- list(dpdbeta0 = dpdbeta0)
  }



  return(derivees_partielles)
}
