#' @title cumLexc_alphaweibull_topred function
#'
#' @description calculates the cumulative excess hazard from a Weibull distribution
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
#' @param theta the parameters of the cumulative excess hazard from a Weibull distribution
#'
#'
#' @param sign_delta only used for mixture cure rate models to specify if the
#' effects or minus the effects of covariates acting on uncured survival to be
#' considered. Default will be sign_delta = "1". The alternative is
#' sign_delta = "-1".
#'
#'
#' @keywords cumLexc_alphaweibull_topred
#'
#' @return This object is a list containing the following components:
#'
#'
#'   \item{cumhaz}{cumulative excess hazard estimates}
#'   \item{usurv}{survival of uncured}
#'   \item{SurvE}{net survival estimates}
#'   \item{cured}{cure fraction}
#'   \item{pt_cure}{the conditionnal probability of being cured knowing they are alive at t}
#'
#' @author Juste Goungounga, Judith Breaud, Eugenie Blandin, Olayidé Boussari, Valérie Jooste
#'
#' @keywords internal

cumLexc_alphaweibull_topred <- function(z_ucured =  z_ucured, z_pcured = z_pcured,
                                        x = x, theta = theta, sign_delta = 1) {
  n_z_pcured <- ncol(z_pcured)
  n_z_ucured <- ncol(z_ucured)
  if (n_z_pcured > 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <-sign_delta*theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]

    pcure <- beta0 + z_pcured %*% betak
    cured <- 1/(1 + exp(-pcure))
    usurv <- (exp(-exp(lambda)*(x)^exp(gamma)))^exp(z_ucured %*% delta)
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1)) * exp(z_ucured %*% delta)
    u_f <- uhaz*usurv
    SurvE <- cured + (1 - cured)*usurv
    cumHazE <- -log(SurvE)
  } else if (n_z_pcured > 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <-sign_delta*theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0 + z_pcured %*% betak
    cured <- 1/(1 + exp(-pcure))
    usurv <- (exp(-exp(lambda)*(x)^exp(gamma)))
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1))
    u_f <- uhaz*usurv
    SurvE <- cured + (1 - cured)*usurv
    cumHazE <- -log(SurvE)

  } else if (n_z_pcured == 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <-sign_delta*theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0
    cured <- as.matrix(rep(1/(1 + exp(-pcure)),length(x)),ncol=1)
    usurv <- (exp(-exp(lambda)*(x)^exp(gamma)))^exp(z_ucured %*% delta)
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1)) * exp(z_ucured %*% delta)
    u_f <- uhaz*usurv
    SurvE <- cured + (1 - cured)*usurv
    cumHazE <- -log(SurvE)
  } else if (n_z_pcured == 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    lambda <- theta[2]
    gamma <- theta[3]
    pcure <- beta0
    cured <- as.matrix(rep(1/(1 + exp(-pcure)),length(x)),ncol=1)
    usurv <- (exp(-exp(lambda)*(x)^exp(gamma)))
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1))
    u_f <- uhaz*usurv
    SurvE <- cured + (1 - cured)*usurv
    cumHazE <- -log(SurvE)
  }
  pt_cure<-cured/SurvE
  return(list(cumhaz = cumHazE,
              cured = cured,
              usurv = usurv,
              SurvE = SurvE,
              pt_cure=pt_cure))
}


