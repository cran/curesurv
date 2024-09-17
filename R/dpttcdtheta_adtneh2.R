#' dpttcdtheta_adtneh2 function
#'
#' @description Partial derivatives of probability to be cure by theta which
#' can be evaluated at t = TTC, from predictions based on non-mixture model
#'  with distribution "tneh".
#'
#' @param object ouput from a model implemented in curesurv
#'
#'
#' @param z_alpha Covariates matrix acting on parameter alpha of the density of
#'  time-to-null excess hazard model
#'
#' @param z_tau Covariates matrix acting on time-to-null parameter.
#'
#' @param x time at which the predictions are provided
#'
#' @param cumLexctopred a pre-prediction parameter obtained with cumLexc_ad2_topred, if NULL will be calculated
#'
#' @param Dpi partial derivative of pi according to theta at time TTC, if NULL will be calculated
#'
#' @param Dsn partial derivative of net survival according to theta at time TTC, if NULL will be calculated
#'
#' @author Juste Goungounga, Judith Breaud, Eugenie Blandin, Olayide Boussari, Valerie Jooste
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

dpttcdtheta_adtneh2 <- function(z_tau,
                               z_alpha,
                               x = x,
                               object,
                               cumLexctopred=NULL,
                               Dpi=NULL,
                               Dsn=NULL) {
  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")
  if(is.null(cumLexctopred)){
    cumLexctopred=cumLexc_ad2_topred(z_tau = z_tau,z_alpha = z_alpha,x=x,theta=object$coefficient)
  }
  if(is.null(Dpi)){
    Dpi <- dpidtheta_adtneh2(z_tau = z_tau,
                                 z_alpha = z_alpha,
                                 x = x,
                                 object,
                                 cumLexctopred=cumLexctopred)
  }
  if(is.null(Dsn)){
    Dsn <- dsndtheta_adtneh2(z_tau = z_tau,
                             z_alpha = z_alpha,
                             x = x,
                             object,
                             cumLexctopred=cumLexctopred,
                             Dpi=Dpi)
  }
  theta <- object$coefficients
  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- theta[1]

  res_pred <- cumLexctopred
  pi <- res_pred$pi
  cumhaz <- res_pred$cumhaz
  netsurv <- res_pred$netsurv
  D <- matrix(0, length(x), length(theta))


  if (n_z_tau == 0 & n_z_alpha == 0) {
    D[, 1] <- (1/netsurv^2) * (Dpi[,1] * netsurv - pi * Dsn[,1])
    D[, 2] <- (1/netsurv^2) * (Dpi[,2] * netsurv - pi * Dsn[,2])
    D[, 3] <- (1/netsurv^2) * (Dpi[,3] * netsurv - pi * Dsn[,3])

  }  else if (n_z_tau > 0 & n_z_alpha > 0) {
    D[, 1] <- (1/netsurv^2) * (Dpi[,1] * netsurv - pi * Dsn[,1])
    D[, 2:(n_z_alpha + 1)] <- (1/matrix(rep(netsurv,n_z_alpha),ncol=n_z_alpha)^2) * (Dpi[,2:(n_z_alpha + 1)] * matrix(rep(netsurv,n_z_alpha),ncol=n_z_alpha) - matrix(rep(pi,n_z_alpha),ncol=n_z_alpha) * Dsn[,2:(n_z_alpha + 1)])


    D[, (n_z_alpha + 2)] <- (1/netsurv^2) * (Dpi[,(n_z_alpha + 2)] * netsurv - pi * Dsn[,(n_z_alpha + 2)])
    D[, (n_z_alpha + 3)] <- (1/netsurv^2) * (Dpi[,(n_z_alpha + 3)] * netsurv - pi * Dsn[,(n_z_alpha + 3)])
    D[, (n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] <- (1/matrix(rep(netsurv,n_z_tau),ncol=n_z_tau)^2) * (Dpi[,(n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] * matrix(rep(netsurv,n_z_tau),ncol=n_z_tau) - matrix(rep(pi,n_z_tau),ncol=n_z_tau) * Dsn[,(n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)])

    }  else if (n_z_tau > 0 & n_z_alpha == 0) {
    D[, 1] <- (1/netsurv^2) * (Dpi[,1] * netsurv - pi * Dsn[,1])
    D[, (n_z_alpha + 2)] <- (1/netsurv^2) * (Dpi[,(n_z_alpha + 2)] * netsurv - pi * Dsn[,(n_z_alpha + 2)])
    D[, (n_z_alpha + 3)] <- (1/netsurv^2) * (Dpi[,(n_z_alpha + 3)] * netsurv - pi * Dsn[,(n_z_alpha + 3)])
    D[, (n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] <- (1/matrix(rep(netsurv,n_z_tau),ncol=n_z_tau)^2) * (Dpi[,(n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] * matrix(rep(netsurv,n_z_tau),ncol=n_z_tau) - matrix(rep(pi,n_z_tau),ncol=n_z_tau) * Dsn[,(n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)])

     }  else if (n_z_tau == 0 & n_z_alpha > 0) {
    D[, 1] <- (1/netsurv^2) * (Dpi[,1] * netsurv - pi * Dsn[,1])
    pi_replicated <- matrix(rep(pi, n_z_alpha), nrow = length(pi), ncol = n_z_alpha)
    D[, 2:(n_z_alpha + 1)] <- (1/netsurv^2) * (Dpi[,2:(n_z_alpha + 1)] * netsurv - pi_replicated * Dsn[,2:(n_z_alpha + 1)])
    D[, (n_z_alpha + 2)] <- (1/netsurv^2) * (Dpi[,(n_z_alpha + 2)] * netsurv - pi * Dsn[,(n_z_alpha + 2)])
    D[, (n_z_alpha + 3)] <- (1/netsurv^2) * (Dpi[,(n_z_alpha + 3)] * netsurv - pi * Dsn[,(n_z_alpha + 3)])
  }

  return(D)
}
