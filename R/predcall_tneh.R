#' @title predcall_tneh function
#'
#' @description calculates the predicted cure indicators from a Time to null
#'  excess hazard model.
#'
#' @param object ouput from a model implemented in curesurv
#'
#' @param pred some predicted estimates
#'
#' @param z_alpha Covariates matrix acting on parameter alpha of the density of
#'  time-to-null excess hazard model
#'
#' @param z_tau Covariates matrix acting on time-to-null parameter.
#'
#' @param x time at which the predictions are provided
#'
#' @param xmax time max at which Pi(t) is calculated.
#'
#' @param level (1-alpha/2)-order quantile of a normal distribution
#'
#' @param epsilon  value fixed by user to estimate the \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
#'   By default \eqn{\epsilon = 0.05}.
#'
#' @param cumLexc_topred  pre prediction obtained either from cumLexc_adtneh2_topred or cumLexc_mul_topred
#'
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


predcall_tneh <- function(object,
                          pred,
                          z_tau = z_tau,
                          z_alpha = z_alpha,
                          x = x,
                          xmax = NULL,
                          level = level,
                          epsilon = epsilon,
                          cumLexc_topred) {
  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")

  if(object$link_tau=="linear"){
    theta <- object$theta

    Dpi <- dpidtheta_adtneh2(z_tau = z_tau,
                             z_alpha = z_alpha,
                             x = x,
                             object,
                             cumLexctopred=cumLexc_topred)


    pi_plain_method <- pi_ic_adtneh2(z_tau = z_tau,
                                     z_alpha = z_alpha,
                                     x = x,
                                     object = object ,
                                     level = level,
                                     cumLexc_topred=cumLexc_topred,
                                     Dpi)
    pred$cured<-pi_plain_method$pi

    pred$lower_bound_pi_tneh <- pi_plain_method$lower_bound

    pred$upper_bound_pi_tneh <- pi_plain_method$upper_bound

    pi_log_method <- pi_ic_adtneh2_log(z_tau = z_tau,
                                       z_alpha = z_alpha,
                                       x = x,
                                       object,
                                       level = level,
                                       cumLexc_topred=cumLexc_topred,
                                       Dpi=Dpi)



    pred$lower_bound_pi_tneh_log <- pi_log_method$lower_bound
    pred$upper_bound_pi_tneh_log <- pi_log_method$upper_bound

    pi_log_log_method <- pi_ic_adtneh2_loglog(z_tau = z_tau,
                                              z_alpha = z_alpha,
                                              x = x,
                                              object,
                                              level = level,
                                              cumLexc_topred=cumLexc_topred,
                                              Dpi=Dpi)

    pred$lower_bound_pi_tneh_loglog <- pi_log_log_method$lower_bound
    pred$upper_bound_pi_tneh_loglog <- pi_log_log_method$upper_bound

    Dexhaz<-dexhazdtheta_adtneh2(z_tau = z_tau,
                                 z_alpha = z_alpha,
                                 x = x,
                                 object,
                                 cumLexctopred=cumLexc_topred)

    exhaz_plain_method<-exhaz_ic_adtneh2(z_tau = z_tau,
                                         z_alpha = z_alpha,
                                         x = x,
                                         object,
                                         level = level,
                                         Dexhaz=Dexhaz)

    pred$lower_bound_exhaz_tneh <- exhaz_plain_method$lower_bound
    pred$upper_bound_exhaz_tneh <- exhaz_plain_method$upper_bound

    exhaz_log_method<-exhaz_ic_adtneh2_log(
      z_tau = z_tau, z_alpha = z_alpha,
      x = x,object,level = level,Dexhaz=Dexhaz)
    pred$lower_bound_exhaz_tneh_log <- exhaz_log_method$lower_bound
    pred$upper_bound_exhaz_tneh_log <- exhaz_log_method$upper_bound

    exhaz_loglog_method<-exhaz_ic_adtneh2_loglog(
      z_tau = z_tau,z_alpha = z_alpha,
      x = x,object,level = level,Dexhaz=Dexhaz)
    pred$lower_bound_exhaz_tneh_loglog <- exhaz_loglog_method$lower_bound
    pred$upper_bound_exhaz_tneh_loglog <- exhaz_loglog_method$upper_bound

    Dsn<-dsndtheta_adtneh2(z_tau=z_tau,
                           z_alpha=z_alpha,
                           x=x,
                           object,
                           cumLexctopred=cumLexc_topred,
                           Dpi=Dpi)

    sn_tneh_plain_method <- sn_ic_adtneh2(z_tau = z_tau,
                                          z_alpha = z_alpha,
                                          x = x,
                                          object,
                                          level = level,
                                          cumLexctopred=cumLexc_topred,
                                          Dsn=Dsn)
    pred$lower_bound_sn_tneh <- sn_tneh_plain_method$lower_bound
    pred$upper_bound_sn_tneh  <- sn_tneh_plain_method$upper_bound

    sn_tneh_log_method <- sn_ic_adtneh2_log(z_tau = z_tau,
                                            z_alpha = z_alpha,
                                            x = x,
                                            object,
                                            level = level,
                                            cumLexctopred=cumLexc_topred,
                                            Dsn=Dsn)
    pred$lower_bound_sn_tneh_log <- sn_tneh_log_method$lower_bound
    pred$upper_bound_sn_tneh_log <- sn_tneh_log_method$upper_bound

    sn_tneh_log_log_method <- sn_ic_adtneh2_loglog(z_tau = z_tau,
                                                   z_alpha = z_alpha,
                                                   x = x,
                                                   object,
                                                   level = level,
                                                   cumLexctopred=cumLexc_topred,
                                                   Dsn=Dsn)
    pred$lower_bound_sn_tneh_loglog <- sn_tneh_log_log_method$lower_bound
    pred$upper_bound_sn_tneh_loglog <- sn_tneh_log_log_method$upper_bound

    Dpt_cure<-dpdtheta_adtneh2(z_alpha = z_alpha,
                                           z_tau = z_tau,
                                           x = x,
                                           object,
                                           cumLexctopred=cumLexc_topred,
                                           Dpi=Dpi,
                                           Dsn=Dsn)
    pt_tneh_plain_method <- pt_cure_ic_adtneh2(z_tau = z_tau,
                                               z_alpha = z_alpha,
                                               x = x,
                                               object,
                                               level = level,
                                               cumLexctopred=cumLexc_topred,
                                               Dpt_cure=Dpt_cure)

    pred$lower_bound_pt_cure_tneh <- pt_tneh_plain_method$lower_bound
    pred$upper_bound_pt_cure_tneh <- pt_tneh_plain_method$upper_bound


    pt_tneh_log_method <- pt_cure_ic_adtneh2_log(z_tau = z_tau,
                                                 z_alpha = z_alpha,
                                                 x = x,
                                                 object,
                                                 level = level,
                                                 cumLexctopred=cumLexc_topred,
                                                 Dpt_cure=Dpt_cure)


    pred$lower_bound_pt_cure_tneh_log <- pt_tneh_log_method$lower_bound
    pred$upper_bound_pt_cure_tneh_log <- pt_tneh_log_method$upper_bound

    pt_tneh_loglog_method <- pt_cure_ic_adtneh2_loglog(z_tau = z_tau,
                                                       z_alpha = z_alpha,
                                                       x = x,
                                                       object,
                                                       level = level,
                                                       cumLexctopred=cumLexc_topred,
                                                       Dpt_cure=Dpt_cure)

    pred$lower_bound_pt_cure_tneh_loglog <- pt_tneh_loglog_method$lower_bound
    pred$upper_bound_pt_cure_tneh_loglog <- pt_tneh_loglog_method$upper_bound


    TTC <-TTC_adtneh2(z_alpha, z_tau, xmax, object, epsilon = epsilon)$TTC
    cumLexc_topred_TTC<-cumLexc_ad2_topred(z_tau,
                                           z_alpha,
                                           x=TTC,
                                           theta=object$coefficient)
    Dpi_TTC<-dpidtheta_adtneh2(z_tau = z_tau,
                               z_alpha = z_alpha,
                               x = TTC,
                               object,
                               cumLexctopred=cumLexc_topred_TTC)
    Dsn_TTC<-dsndtheta_adtneh2(z_tau = z_tau,
                               z_alpha = z_alpha,
                               x = TTC,
                               object,
                               cumLexctopred=cumLexc_topred_TTC,
                               Dpi = Dpi_TTC)

    Dpttc<-dpttcdtheta_adtneh2(z_tau,
                                    z_alpha,
                                    x = x,
                                    object,
                                    cumLexctopred=cumLexc_topred_TTC,
                                    Dpi=Dpi_TTC,
                                    Dsn=Dsn_TTC)

    Var_TTC<-var_TTC_tneh2(z_alpha, z_tau, xmax, object, epsilon = epsilon,TTC,
                                     cumLexc_topred_TTC,
                                     Dpttc)


    ttc_ic_tneh_plain_method <- TTC_ic_adtneh2(z_alpha, z_tau, xmax, object,
                                               epsilon = epsilon, level = level,
                                               TTC = TTC,
                                               varTTC = Var_TTC)

    pred$time_to_cure_ttc <- TTC

    pred$lower_bound_TTC_tneh <- ttc_ic_tneh_plain_method$lower_bound
    pred$upper_bound_TTC_tneh  <- ttc_ic_tneh_plain_method$upper_bound

    ttc_ic_tneh_log_method <- TTC_ic_adtneh2_log(z_alpha, z_tau, xmax, object,
                                                 epsilon = epsilon, level = level,
                                                 TTC = TTC,varTTC = Var_TTC)
    pred$lower_bound_TTC_tneh_log <- ttc_ic_tneh_log_method$lower_bound
    pred$upper_bound_TTC_tneh_log <- ttc_ic_tneh_log_method$upper_bound

    pred$model <- object$model
    pred$dist <- object$dist
    pred$link_tau <- object$link_tau
    pred$epsilon <- epsilon


    return(pred)

  }else if(object$link_tau=="loglinear"){
    theta <- object$theta


    Dpi<-dpidtheta_multneh(z_tau = z_tau,
                                       z_alpha = z_alpha,
                                       x = x,
                                       object,
                                       cumLexctopred=cumLexc_topred)
    pi_plain_method <- pi_ic_multneh(z_tau = z_tau,
                                     z_alpha = z_alpha,
                                     x = x,
                                     object = object ,
                                     level = level,
                                     cumLexctopred = cumLexc_topred,
                                     Dpi=Dpi)
    pred$cured<-pi_plain_method$pi
    pred$lower_bound_pi_tneh <- pi_plain_method$lower_bound
    pred$upper_bound_pi_tneh <- pi_plain_method$upper_bound

    pi_log_method <- pi_ic_multneh_log(z_tau = z_tau,
                                       z_alpha = z_alpha,
                                       x = x,
                                       object,
                                       level = level,cumLexctopred = cumLexc_topred,
                                       Dpi=Dpi)

    pred$lower_bound_pi_tneh_log <- pi_log_method$lower_bound
    pred$upper_bound_pi_tneh_log <- pi_log_method$upper_bound

    pi_log_log_method <- pi_ic_multneh_loglog(z_tau = z_tau,
                                              z_alpha = z_alpha,
                                              x = x,
                                              object,
                                              level = level,
                                              cumLexctopred = cumLexc_topred,
                                              Dpi=Dpi)

    pred$lower_bound_pi_tneh_loglog <- pi_log_log_method$lower_bound
    pred$upper_bound_pi_tneh_loglog <- pi_log_log_method$upper_bound


    Dexhaz<-dexhazdtheta_multneh(z_tau = z_tau,
                                 z_alpha = z_alpha,
                                 x = x,
                                 object,
                                 cumLexctopred=cumLexc_topred)

    exhaz_plain_method<-exhaz_ic_multneh(z_tau = z_tau,
                                         z_alpha = z_alpha,
                                         x = x,
                                         object,
                                         level = level,
                                         Dexhaz=Dexhaz)

    pred$lower_bound_exhaz_tneh <- exhaz_plain_method$lower_bound
    pred$upper_bound_exhaz_tneh <- exhaz_plain_method$upper_bound

    exhaz_log_method<-exhaz_ic_multneh_log(
      z_tau = z_tau, z_alpha = z_alpha,
      x = x,object,level = level,Dexhaz=Dexhaz)
    pred$lower_bound_exhaz_tneh_log <- exhaz_log_method$lower_bound
    pred$upper_bound_exhaz_tneh_log <- exhaz_log_method$upper_bound

    exhaz_loglog_method<-exhaz_ic_multneh_loglog(
      z_tau = z_tau,z_alpha = z_alpha,
      x = x,object,level = level,Dexhaz=Dexhaz)
    pred$lower_bound_exhaz_tneh_loglog <- exhaz_loglog_method$lower_bound
    pred$upper_bound_exhaz_tneh_loglog <- exhaz_loglog_method$upper_bound


    Dsn<-dsndtheta_multneh(z_tau, z_alpha, x, object,cumLexctopred=cumLexc_topred,Dpi)
    sn_tneh_plain_method <- sn_ic_multneh(z_tau = z_tau,
                                          z_alpha = z_alpha,
                                          x = x,
                                          object,
                                          level = level,
                                          cumLexctopred = cumLexc_topred,
                                          Dsn=Dsn)
    pred$lower_bound_sn_tneh <- sn_tneh_plain_method$lower_bound
    pred$upper_bound_sn_tneh  <- sn_tneh_plain_method$upper_bound

    sn_tneh_log_method <- sn_ic_multneh_log(z_tau = z_tau,
                                            z_alpha = z_alpha,
                                            x = x,
                                            object,
                                            level = level,
                                            cumLexctopred = cumLexc_topred,
                                            Dsn=Dsn)
    pred$lower_bound_sn_tneh_log <- sn_tneh_log_method$lower_bound
    pred$upper_bound_sn_tneh_log <- sn_tneh_log_method$upper_bound

    sn_tneh_log_log_method <- sn_ic_multneh_loglog(z_tau = z_tau,
                                                   z_alpha = z_alpha,
                                                   x = x,
                                                   object,
                                                   level = level,
                                                   cumLexctopred = cumLexc_topred,
                                                   Dsn=Dsn)
    pred$lower_bound_sn_tneh_loglog <- sn_tneh_log_log_method$lower_bound
    pred$upper_bound_sn_tneh_loglog <- sn_tneh_log_log_method$upper_bound

    Dpt_cure<-dptdtheta_multneh(z_alpha = z_alpha,
                                          z_tau = z_tau,
                                          x = x,
                                          object,
                                          cumLexctopred=cumLexc_topred,
                                          Dpi=Dpi,
                                          Dsn=Dsn)

    pt_tneh_plain_method <- pt_cure_ic_multneh(z_tau = z_tau,
                                               z_alpha = z_alpha,
                                               x = x,
                                               object,
                                               level = level,
                                               cumLexctopred = cumLexc_topred,
                                               Dpt_cure=Dpt_cure)


    pred$lower_bound_pt_cure_tneh <- pt_tneh_plain_method$lower_bound
    pred$upper_bound_pt_cure_tneh <- pt_tneh_plain_method$upper_bound
    pt_tneh_log_method <- pt_cure_ic_multneh_log(z_tau = z_tau,
                                                 z_alpha = z_alpha,
                                                 x = x,
                                                 object,
                                                 level = level,
                                                 cumLexctopred = cumLexc_topred,
                                                 Dpt_cure=Dpt_cure)


    pred$lower_bound_pt_cure_tneh_log <- pt_tneh_log_method$lower_bound
    pred$upper_bound_pt_cure_tneh_log <- pt_tneh_log_method$upper_bound

    pt_tneh_loglog_method <- pt_cure_ic_multneh_loglog(z_tau = z_tau,
                                                       z_alpha = z_alpha,
                                                       x = x,
                                                       object,
                                                       level = level,
                                                       cumLexctopred = cumLexc_topred,
                                                       Dpt_cure=Dpt_cure)

    pred$lower_bound_pt_cure_tneh_loglog <- pt_tneh_loglog_method$lower_bound
    pred$upper_bound_pt_cure_tneh_loglog <- pt_tneh_loglog_method$upper_bound

    TTC <-TTC_multneh(z_alpha, z_tau, xmax, object, epsilon = epsilon)$TTC
    cumLexctopred_TTC<-cumLexc_mul_topred(z_tau,
                                      z_alpha,
                                      x = TTC,
                                      object$coefficient)
    Dpi_TTC<- dpidtheta_multneh(z_tau = z_tau,
                                z_alpha = z_alpha,
                                x = TTC,object,cumLexctopred = cumLexctopred_TTC)
    Dsn_TTC<- dsndtheta_multneh(z_tau = z_tau,
                                z_alpha = z_alpha,
                                x = TTC,object,cumLexctopred = cumLexctopred_TTC,Dpi = Dpi_TTC)
    DpTTC<-dpttcdtheta_multneh(z_tau,z_alpha,x = x, object,
                               res_pred=cumLexctopred_TTC,Dpi=Dpi_TTC,Dsn=Dsn_TTC)


    ttc_ic_tneh_plain_method <- TTC_ic_multneh(z_alpha, z_tau, xmax, object,
                                               epsilon = epsilon, level = level,
                                               TTC=TTC,
                                               DpTTC,
                                               cumLexctopred=cumLexctopred_TTC)

    pred$time_to_cure_ttc <- TTC
    pred$lower_bound_TTC_tneh <- ttc_ic_tneh_plain_method$lower_bound
    pred$upper_bound_TTC_tneh  <- ttc_ic_tneh_plain_method$upper_bound

    ttc_ic_tneh_log_method <- TTC_ic_multneh_log(z_alpha, z_tau, xmax, object,
                                                 epsilon = epsilon, level = level,
                                                 TTC=TTC,
                                                 DpTTC,
                                                 cumLexctopred=cumLexctopred_TTC)
    pred$lower_bound_TTC_tneh_log <- ttc_ic_tneh_log_method$lower_bound
    pred$upper_bound_TTC_tneh_log <- ttc_ic_tneh_log_method$upper_bound

    pred$model <- object$model
    pred$dist <- object$dist
    pred$link_tau <- object$link_tau
    pred$epsilon <- epsilon
    return(pred)
  }else{
    stop("Not implemented yet for this link")
  }



  }
