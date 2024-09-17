#' @title Prediction for a curesurv cure model
#'
#' @description return predicted (excess) hazard, (net) survival, cure fraction and
#' time to null excess hazard or time to cure.
#'
#' @param object Output from \code{curesurv} function
#'
#'
#' @param newdata the new data to be specified for predictions; If else,
#'  predictions are made using the data provided during the estimation step in
#'  order to obtain the output from curesurv function.
#'
#' @param xmax maximum time at which Time-to-Cure is evaluated numerically.
#'
#' @param level \eqn{1-\frac{\alpha}{2}}-order quantile of a normal distribution for the confidence intervals
#'
#' @param epsilon  value fixed by user to estimate the TTC \eqn{Pi(t) \ge 1-\epsilon}.
#'  By default epsilon = 0.05.
#'
#' @param sign_delta sign of effect of delta on covariates acting on survival function,
#'  positive by default "sign_delta = 1" and alternative is "sign_delta = -1"
#'
#'
#' @param ... additional parameters
#'
#' @keywords predict.curesurv
#'
#' @return An object of class \code{c("pred_curesurv", "data.frame")}.
#' This object is a list containing the following components:
#'
#'
#'   \item{time}{time in the input new data}
#'
#'   \item{ex_haz}{predicted excess hazard at the time provided in the new data}
#'
#'   \item{netsurv}{predicted net survival at the time provided in the new data}
#'
#'   \item{pt_cure}{probability to be cured}
#'
#'   \item{tau}{time to null in model TNEH when object corresponds to the
#'   results from Boussari model or its extension.}
#'
#'
#'   \item{netsurv_tau}{pi or net survival at time tau when object corresponds to the
#'   results from Boussari model or its extension.}
#'
#'   \item{time_to_cure_ttc}{time to cure (TTC)}
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
#' @examples
#'
#' library("curesurv")
#' library("survival")
#'
#' fit_m2_ml <- curesurv(Surv(time_obs, event) ~ age_cr|age_cr,
#'                    pophaz = "ehazard",
#'                    cumpophaz = "cumehazard",
#'                    model = "mixture",
#'                    data = pancreas_data,
#'                    method_opt = "L-BFGS-B")
#'
#'  fit_m2_ml
#'
#'  newdata <- pancreas_data[2,]
#'
#'  predict(object = fit_m2_ml, newdata = newdata)
#'
#' ## Non mixture cure model
#' ### TNEH model
#'
#' #### Additive parametrization
#'
#' testiscancer$age_crmin <- (testiscancer$age- min(testiscancer$age)) /
#'               sd(testiscancer$age)
#'
#' fit_m1_ad_tneh <- curesurv(Surv(time_obs, event) ~ z_tau(age_crmin) +
#'                           z_alpha(age_crmin),
#'                           pophaz = "ehazard",
#'                           cumpophaz = "cumehazard",
#'                           model = "nmixture", dist = "tneh",
#'                           link_tau = "linear",
#'                           data = testiscancer,
#'                           method_opt = "L-BFGS-B")
#'
#'  fit_m1_ad_tneh
#'
#'  predict(object = fit_m1_ad_tneh, newdata = testiscancer[3:6,])
#'
#'  #mean of age
#'  newdata1 <- with(testiscancer,
#'  expand.grid(event = 0, age_crmin = mean(age_crmin), time_obs  = seq(0.001,10,0.1)))
#'
#'  pred_agemean <- predict(object = fit_m1_ad_tneh, newdata = newdata1)
#'
#'
#'  #max of age
#'  newdata2 <- with(testiscancer,
#'  expand.grid(event = 0,
#'  age_crmin = max(age_crmin),
#'   time_obs  = seq(0.001,10,0.1)))
#'
#'  pred_agemax <- predict(object = fit_m1_ad_tneh, newdata = newdata2)
#'  head(pred_agemax)
#'
#'
#'
#
#'
#' @seealso [curesurv::print.curesurv()], [curesurv::curesurv()], `browseVignettes("curesurv")`
#'
#' @import survival
#'
#' @export
predict.curesurv <- function(object,
                             newdata = NULL,
                             xmax = 10^9,
                             level = 0.975,
                             epsilon = 0.05, sign_delta = 1, ...) {

  pred <- predictform(object, newdata,
                      xmax,
                      level,
                      epsilon,
                      sign_delta = 1,
                      ...)

  return(pred)
}

