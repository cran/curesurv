#' @title dexhazdtheta_wei function
#'
#' @description Produce partial derivatives of excess hazard
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
#' @param cumLexctopred pre prediction obtained with cumLexc_alphaweibull_topred
#'
#' @keywords internal

dexhazdtheta_wei <- function(z_pcured = z_pcured,
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

    dehdbeta0 <- (exp(-pcure)/1+exp(-pcure))*(-cured*u_f*(cured+(1-cured))*usurv-(1-usurv)*(1-cured)*u_f)/(SurvE^2)

    dehdbeta_star <- matrix(rep(dehdbeta0,n_z_pcured),ncol=n_z_pcured) * z_pcured

    dsudlambda<--exp(lambda)*exp(z_ucured%*%delta)*x^(exp(gamma))*usurv
    dfudlambda<-u_f*(1-exp(lambda)*exp(z_ucured%*%delta)*x^(exp(gamma)))
    dehdlambda <- ((1-cured)*dfudlambda*SurvE-(1-cured)^2*dsudlambda*u_f)/(SurvE^2)

    dsudgamma<- -exp(lambda)*exp(z_ucured%*%delta)*log(x)*exp(gamma)*x^(exp(gamma))*usurv
    dfudgamma<-u_f*(1+exp(gamma)*log(x))+dsudgamma*u_f/usurv
    dehdgamma <- ((1-cured)*dfudgamma*SurvE-(1-cured)^2*dsudgamma*u_f)/(SurvE^2)

    dehddelta <- z_ucured * matrix(rep(dehdgamma,n_z_ucured),ncol=n_z_ucured)


    derivees_partielles <- list(dehdbeta0 = dehdbeta0,
                                dehdbeta_star = dehdbeta_star,
                                dehdlambda = dehdlambda,
                                dehdgamma = dehdgamma,
                                dehddelta = dehddelta)

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

    dehdbeta0 <- (exp(-pcure)/1+exp(-pcure))*(-cured*u_f*(cured+(1-cured))*usurv-(1-usurv)*(1-cured)*u_f)/(SurvE^2)

    dehdbeta_star <- matrix(rep(dehdbeta0,n_z_pcured),ncol=n_z_pcured) * z_pcured

    dsudlambda<--exp(lambda)*1*x^(exp(gamma))*usurv
    dfudlambda<-u_f*(1-exp(lambda)*1*x^(exp(gamma)))
    dehdlambda <- ((1-cured)*dfudlambda*SurvE-(1-cured)^2*dsudlambda*u_f)/(SurvE^2)

    dsudgamma<- -exp(lambda)*1*log(x)*exp(gamma)*x^(exp(gamma))*usurv
    dfudgamma<-u_f*(1+exp(gamma)*log(x))+dsudgamma*u_f/usurv
    dehdgamma <- ((1-cured)*dfudgamma*SurvE-(1-cured)^2*dsudgamma*u_f)/(SurvE^2)

    derivees_partielles <- list(dehdbeta0 = dehdbeta0,
                                dehdbeta_star = dehdbeta_star,
                                dehdlambda = dehdlambda,
                                dehdgamma = dehdgamma)

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

    dehdbeta0 <- (exp(-pcure)/1+exp(-pcure))*(-cured*u_f*(cured+(1-cured))*usurv-(1-usurv)*(1-cured)*u_f)/(SurvE^2)

    dsudlambda<--exp(lambda)*exp(z_ucured%*%delta)*x^(exp(gamma))*usurv
    dfudlambda<-u_f*(1-exp(lambda)*exp(z_ucured%*%delta)*x^(exp(gamma)))
    dehdlambda <- ((1-cured)*dfudlambda*SurvE-(1-cured)^2*dsudlambda*u_f)/(SurvE^2)

    dsudgamma<- -exp(lambda)*exp(z_ucured%*%delta)*log(x)*exp(gamma)*x^(exp(gamma))*usurv
    dfudgamma<-u_f*(1+exp(gamma)*log(x))+dsudgamma*u_f/usurv
    dehdgamma <- ((1-cured)*dfudgamma*SurvE-(1-cured)^2*dsudgamma*u_f)/(SurvE^2)

    dehddelta <- z_ucured * matrix(rep(dehdgamma,n_z_ucured),ncol=n_z_ucured)

    derivees_partielles <- list(dehdbeta0 = dehdbeta0,
                                dehdlambda = dehdlambda,
                                dehdgamma = dehdgamma,
                                dehddelta = dehddelta)

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

    dehdbeta0 <- (exp(-pcure)/1+exp(-pcure))*(-cured*u_f*(cured+(1-cured))*usurv-(1-usurv)*(1-cured)*u_f)/(SurvE^2)

    dsudlambda<--exp(lambda)*1*x^(exp(gamma))*usurv
    dfudlambda<-u_f*(1-exp(lambda)*1*x^(exp(gamma)))
    dehdlambda <- ((1-cured)*dfudlambda*SurvE-(1-cured)^2*dsudlambda*u_f)/(SurvE^2)

    dsudgamma<- -exp(lambda)*1*log(x)*exp(gamma)*x^(exp(gamma))*usurv
    dfudgamma<-u_f*(1+exp(gamma)*log(x))+dsudgamma*u_f/usurv
    dehdgamma <- ((1-cured)*dfudgamma*SurvE-(1-cured)^2*dsudgamma*u_f)/(SurvE^2)

    derivees_partielles <- list(dehdbeta0 = dehdbeta0,
                                dehdlambda = dehdlambda,
                                dehdgamma = dehdgamma)
  }

  return(derivees_partielles)
}
