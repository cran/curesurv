#' @title var_TTC_tneh2 function
#'
#' @description calculates the variance of TTC in a non-mixture model with
#' distribution "TNEH"
#'
#'
#' @param object ouput from a model implemented in curesurv
#'
#'
#' @param z_alpha Covariates matrix acting on parameter alpha of the density of
#'  time-to-null excess hazard model
#'
#' @param z_tau Covariates matrix acting on time-to-null parameter.
#'
#' @param xmax time max at which Pi(t) is calculated.
#'
#' @param epsilon value fixed by user to estimate the TTC \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
#'   By default \eqn{\epsilon = 0.05}.
#'
#' @param TTC time to cure calculated by TTC_adtneh2
#'
#' @param cumLexctopred pre prediction, if NULL recalculated
#'
#' @param Dpttc Partial derivatives of probability to be cure by theta which
#' can be evaluated at t = TTC, if NULL it is recalculated
#'
#' @author Juste Goungounga, Judith Breaud, Eugenie Blandin, Olayidé Boussari, Valérie Jooste
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

var_TTC_tneh2 <- function(z_alpha, z_tau, xmax, object, epsilon = epsilon,TTC,
                          cumLexctopred=NULL,Dpttc=NULL) {


  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")


  theta <- object$coefficients
  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- theta[1]

  if (n_z_tau == 0 & n_z_alpha == 0) {

    alpha <- theta[1]
    beta <- theta[2]
    tau <- theta[3]

  }


  else if (n_z_tau > 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- theta[n_z_alpha + 2]
    beta2 <- beta
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    alpha <- alpha0 + z_alpha %*% alpha_k
    tau <- tau0 + z_tau %*% tau_z

  }


  else if (n_z_tau > 0 & n_z_alpha == 0) {

    alpha <- theta[1]
    beta <- theta[n_z_alpha + 2]
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    tau <- tau0 + z_tau %*% tau_z

  }


  else if (n_z_tau == 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- theta[n_z_alpha + 2]
    tau0 <- theta[n_z_alpha + 2 + 1]
    alpha <- alpha0 + z_alpha %*% alpha_k
    tau <- tau0

  }
  if(is.null(cumLexctopred)){
    cumLexctopred <- cumLexc_ad2_topred(z_tau,
                                        z_alpha,
                                        x = TTC,
                                        theta)
  }

  pt_cure <- cumLexctopred$pt_cure

  if(is.null(Dpttc)){
    Dpttc <- dpttcdtheta_adtneh2(z_tau,
                                       z_alpha, x = TTC, object,cumLexctopred)

  }



  var_pt_cure <- Dpttc %*% object$varcov_star %*% t(Dpttc)


  dpdt <- ((TTC/tau)^(alpha - 1)) * ((1 - (TTC/tau))^(beta - 1)) * pt_cure
  var_ttc <- dpdt^(-2) * diag(var_pt_cure)


  return(var_ttc)

}
