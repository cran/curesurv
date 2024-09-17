#' @title predcall_wei function
#'
#' @description calculates the predicted cure indicators from a mixture cure
#'  model with the survival of uncured specified by a Weibull distribution.
#'
#'
#' @param object ouput from a model implemented using curesurv
#'
#' @param pred some predicted estimates
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param x time at which the predictions are provided
#'
#' @param level (1-alpha/2)-order quantile of a normal distribution
#'
#' @param epsilon  value fixed by user to estimate the TTC \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
#'   By default \eqn{\epsilon = 0.05}.
#'
#' @param  sign_delta only used for mixture cure rate models to specify if the
#' effects or minus the effects of covariates acting on uncured survival to be
#' considered. Default will be sign_delta = "1". The alternative is
#' sign_delta = "-1".
#'
#' @param cumLexctopred pre prediction obtained by cumLexc_alphaweibull_topred
#'
#'
#' @author Juste Goungounga, Judith Breaud, Olayide Boussari, Valerie Jooste
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
#'
#'  Phillips N, Coldman A, McBride ML. Estimating cancer prevalence using
#'  mixture models for cancer survival. Stat Med. 2002 May 15;21(9):1257-70.
#'  doi: 10.1002/sim.1101. PMID: 12111877.
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/12111877/}{pubmed})
#'
#'
#'   De Angelis R, Capocaccia R, Hakulinen T, Soderman B, Verdecchia A. Mixture
#'   models for cancer survival analysis: application to population-based data
#'   with covariates. Stat Med. 1999 Feb 28;18(4):441-54.
#'   doi: 10.1002/(sici)1097-0258(19990228)18:4<441::aid-sim23>3.0.co;2-m.
#'   PMID: 10070685.
#'   (\href{https://pubmed.ncbi.nlm.nih.gov/10070685/}{pubmed})
#'
#' @keywords internal

predcall_wei <- function(object,
                         pred,
                         z_pcured = z_pcured,
                         z_ucured = z_ucured,
                         x = x,
                         level = level,
                         epsilon = epsilon, sign_delta = 1,
                         cumLexctopred) {

  theta <- object$coefficients

  pred$time_to_cure_ttc <- TTC_wei(z_pcured = z_pcured,
                                   z_ucured = z_ucured,
                                   theta = object$coefficients,
                                   epsilon = epsilon)

  pred$p_ttc <- pt_cure_wei(z_pcured = z_pcured, z_ucured = z_ucured,
                            x = pred$time_to_cure_ttc, theta = object$coefficients)


  cumLexc_topred_TTC<-cumLexc_alphaweibull_topred(z_ucured =  z_ucured,
                              z_pcured = z_pcured,
                              x = pred$time_to_cure_ttc,
                              theta,  sign_delta = sign_delta)

  pred$netsurv_ttc <- exp(-cumLexc_topred_TTC$cumhaz)

  Dpt_cure <-dpdtheta_wei(
        z_pcured = z_pcured,
        z_ucured = z_ucured,
        x = x,
        theta = object$coefficients,
        cumLexctopred
      )

  pt_ic <- pt_cure_ic_wei(object,
                          z_pcured = z_pcured,
                          z_ucured = z_ucured,
                          x = x,
                          level = level,cumLexctopred,Dpt_cure)


  pred$lower_bound_pt_ic <- pt_ic$lower_bound
  pred$upper_bound_pt_ic <- pt_ic$upper_bound

  Dpt_cure_TTC <-dpdtheta_wei(
    z_pcured = z_pcured,
    z_ucured = z_ucured,
    x = pred$time_to_cure,
    theta = object$coefficients,
    cumLexc_topred_TTC
  )

  pttc_ic <- pt_cure_ic_wei(object, z_pcured = z_pcured,
                            z_ucured = z_ucured,
                            x = pred$time_to_cure,
                            level,cumLexc_topred_TTC,Dpt_cure_TTC)
  pred$lower_bound_pttc_ic <- pttc_ic$lower_bound
  pred$upper_bound_pttc_ic <- pttc_ic$upper_bound

  pt_ic_log <- pt_cure_ic_log_wei(object, z_pcured = z_pcured,
                                  z_ucured = z_ucured,
                                  x= x,
                                  level,cumLexctopred,Dpt_cure)

  pred$lower_bound_pt_ic_log <- pt_ic_log$lower_bound
  pred$upper_bound_pt_ic_log <- pt_ic_log$upper_bound

  pttc_ic_log <- pt_cure_ic_log_wei(object, z_pcured = z_pcured,
                                    z_ucured = z_ucured,
                                    x = pred$time_to_cure,
                                    level,cumLexc_topred_TTC,Dpt_cure_TTC)
  pred$lower_bound_pttc_ic_log <- pttc_ic_log$lower_bound
  pred$upper_bound_pttc_ic_log <- pttc_ic_log$upper_bound
  pt_ic_loglog <- pt_cure_ic_loglog_wei(object, z_pcured = z_pcured,
                                        z_ucured = z_ucured,
                                        x = x,
                                        level,
                                        cumLexctopred,Dpt_cure)

  pred$lower_bound_pt_ic_loglog <- pt_ic_loglog$lower_bound
  pred$upper_bound_pt_ic_loglog <- pt_ic_loglog$upper_bound

  pttc_ic_loglog <- pt_cure_ic_loglog_wei(object, z_pcured = z_pcured,
                                          z_ucured = z_ucured,
                                          x = pred$time_to_cure,
                                          level, cumLexc_topred_TTC,Dpt_cure_TTC)
  pred$lower_bound_pttc_ic_loglog <- pttc_ic_loglog$lower_bound
  pred$upper_bound_pttc_ic_loglog  <- pttc_ic_loglog$upper_bound

  Dexhaz<-dexhazdtheta_wei(z_pcured = z_pcured,
                           z_ucured = z_ucured,
                           x = x, theta,cumLexctopred)
  exhaz_ic <- exhaz_ic_wei(object,
                     z_pcured = z_pcured,
                     z_ucured = z_ucured,
                     x = x,
                     level,cumLexctopred,Dexhaz,exhaz=pred$ex_haz)

  pred$lower_bound_exhaz_ic <- exhaz_ic$lower_bound
  pred$upper_bound_exhaz_ic <- exhaz_ic$upper_bound

  exhaz_ic_log <- exhaz_ic_wei_log(object,
                           z_pcured = z_pcured,
                           z_ucured = z_ucured,
                           x = x,
                           level,cumLexctopred,Dexhaz,exhaz=pred$ex_haz)

  pred$lower_bound_exhaz_ic_log <- exhaz_ic_log$lower_bound
  pred$upper_bound_exhaz_ic_log <- exhaz_ic_log$upper_bound
  exhaz_ic_loglog <- exhaz_ic_wei_loglog(object,
                                   z_pcured = z_pcured,
                                   z_ucured = z_ucured,
                                   x = x,
                                   level,cumLexctopred,Dexhaz,exhaz=pred$ex_haz)

  pred$lower_bound_exhaz_ic_loglog <- exhaz_ic_loglog$lower_bound
  pred$upper_bound_exhaz_ic_loglog <- exhaz_ic_loglog$upper_bound

  Dpi<-dpidtheta_wei(
    z_ucured =  z_ucured,
    z_pcured = z_pcured,
    x = x,
    theta = object$coefficients,
    cumLexctopred
  )
  pi_ic <- pi_ic_wei(object,
                     z_pcured = z_pcured,
                     z_ucured = z_ucured,
                     x = x,
                     level,cumLexctopred,Dpi)

  pred$lower_bound_pi_ic <- pi_ic$lower_bound
  pred$upper_bound_pi_ic <- pi_ic$upper_bound

  pi_ic_log <- pi_ic_log_wei(object,
                             z_pcured = z_pcured,
                             z_ucured = z_ucured,
                             x = x,
                             level,cumLexctopred,Dpi)

  pred$lower_bound_pi_ic_log <- t(pi_ic_log$lower_bound)
  pred$upper_bound_pi_ic_log <- t(pi_ic_log$upper_bound)
  pi_ic_loglog <- pi_ic_loglog_wei(object,
                             z_pcured = z_pcured,
                             z_ucured = z_ucured,
                             x = x,
                             level,cumLexctopred,Dpi)

  pred$lower_bound_pi_ic_loglog <- pi_ic_loglog$lower_bound
  pred$upper_bound_pi_ic_loglog <- pi_ic_loglog$upper_bound
  Dsn<-dsndtheta_wei(
    z_ucured =  z_ucured,
    z_pcured = z_pcured,
    x = x,
    theta = object$coefficients,
    cumLexctopred
  )
  sn_ic <- sn_ic_wei(object,
                     z_pcured = z_pcured,
                     z_ucured = z_ucured,
                     x = x,
                     level,cumLexctopred,Dsn)

  pred$lower_bound_sn_ic <- sn_ic$lower_bound
  pred$upper_bound_sn_ic <- sn_ic$upper_bound

  Dsn_TTC<-dsndtheta_wei(
    z_ucured =  z_ucured,
    z_pcured = z_pcured,
    x = pred$time_to_cure_ttc,
    theta = object$coefficients,
    cumLexc_topred_TTC
  )
  snttc_ic <- sn_ic_wei(object,
                        z_pcured = z_pcured,
                        z_ucured = z_ucured,
                        x = pred$time_to_cure,
                        level,cumLexc_topred_TTC,Dsn_TTC)

  pred$lower_bound_snttc_ic <- snttc_ic$lower_bound
  pred$upper_bound_snttc_ic <- snttc_ic$upper_bound

  sn_ic_log <- sn_ic_log_wei(object,
                         z_pcured = z_pcured,
                         z_ucured = z_ucured,
                         x = x,
                         level,cumLexctopred,Dsn)

  pred$lower_bound_sn_ic_log <- sn_ic_log$lower_bound
  pred$upper_bound_sn_ic_log <- sn_ic_log$upper_bound

  snttc_ic_log <- sn_ic_log_wei(object,
                                z_pcured = z_pcured,
                                z_ucured = z_ucured,
                                x = pred$time_to_cure,
                                level,
                                cumLexc_topred_TTC,Dsn_TTC)

  pred$lower_bound_snttc_ic_log <- snttc_ic_log$lower_bound
  pred$upper_bound_snttc_ic_log <- snttc_ic_log$upper_bound

  logCumHaz_ic <- logCumHaz_ic_wei(object,
                               z_pcured = z_pcured,
                               z_ucured = z_ucured,
                               x = x,
                               level,cumLexctopred)

  pred$lower_bound_logCumHaz_ic <- logCumHaz_ic$lower_bound
  pred$upper_bound_logCumHaz_ic <- logCumHaz_ic$upper_bound
  pred$lower_bound_sn_ic_loglog <- exp(-exp(logCumHaz_ic$upper_bound))
  pred$upper_bound_sn_ic_loglog <- exp(-exp(logCumHaz_ic$lower_bound))

  logCumHazttc_ic <- logCumHaz_ic_wei(object,
                                  z_pcured = z_pcured,
                                  z_ucured = z_ucured,
                                  x = pred$time_to_cure,
                                  level,cumLexc_topred_TTC)

  pred$lower_bound_snttc_ic_loglog <- exp(-exp(logCumHazttc_ic$upper_bound))
  pred$upper_bound_snttc_ic_loglog <- exp(-exp(logCumHazttc_ic$lower_bound))

  Dttc<-dTTCdtheta_wei(z_ucured =  z_ucured,
                       z_pcured = z_pcured,
                       theta = object$coefficients,
                       epsilon = epsilon, TTC=pred$time_to_cure_ttc)

  TTC_ic <- TTC_ic_wei(object, z_pcured = z_pcured,
                       z_ucured = z_ucured,
                       epsilon = epsilon, level,TTC=pred$time_to_cure_ttc,Dttc)

  pred$lower_bound_TTC_ic <- TTC_ic$lower_bound
  pred$upper_bound_TTC_ic <- TTC_ic$upper_bound

  TTC_ic_log <- TTC_ic_log_wei(object, z_pcured = z_pcured,
                               z_ucured = z_ucured,
                               epsilon = epsilon, level,TTC=pred$time_to_cure_ttc,Dttc)
  pred$lower_bound_TTC_ic_log <- TTC_ic_log$lower_bound
  pred$upper_bound_TTC_ic_log  <- TTC_ic_log$upper_bound


  if(ncol(z_pcured)==0&ncol(z_ucured)==0){
    pred$var_TTC_Jakobsen <- var_TTC_Jakobsen_wei(object,
                                                  z_ucured =  z_ucured,
                                                  z_pcured = z_pcured,
                                                  epsilon = epsilon,
                                                  TTC=pred$time_to_cure_ttc,
                                                  cumLexc_topred_TTC)

  }else{
    pred$var_TTC_Jakobsen <- diag(var_TTC_Jakobsen_wei(object,
                                                       z_ucured =  z_ucured,
                                                       z_pcured = z_pcured,
                                                       epsilon = epsilon,
                                                       TTC=pred$time_to_cure_ttc,
                                                       cumLexc_topred_TTC))

  }

  pred$var_TTC <- diag(var_TTC_wei(object,
                         z_ucured =  z_ucured,
                         z_pcured = z_pcured,
                         epsilon = epsilon,TTC=pred$time_to_cure_ttc))


  pred$varlogTTC <- diag(varlogTTC_wei(object,
                                       z_ucured =  z_ucured,
                                       z_pcured = z_pcured,
                                       epsilon = epsilon,TTC=pred$time_to_cure_ttc))
  if(ncol(z_pcured)==0&ncol(z_ucured)==0){
    pred$varlogTTC_Jakobsen <- varlogTTC_Jakobsen_wei(object,
                                                           z_ucured =  z_ucured,
                                                           z_pcured = z_pcured,
                                                           epsilon = epsilon,
                                                      TTC=pred$time_to_cure_ttc,
                                                      cumLexc_topred_TTC)

  }else{
    pred$varlogTTC_Jakobsen <- diag(varlogTTC_Jakobsen_wei(object,
                                                           z_ucured =  z_ucured,
                                                           z_pcured = z_pcured,
                                                           epsilon = epsilon,
                                                           TTC=pred$time_to_cure_ttc,
                                                           cumLexc_topred_TTC))

  }


  pred$model <- object$model
  pred$dist <- object$dist
  pred$epsilon <- epsilon

  return(pred)
}
