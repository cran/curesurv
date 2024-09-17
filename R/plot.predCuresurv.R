#' @title plot method for curesurv prediction objects
#'
#' @description Produces figures of (excess) hazard, (net) survival and
#' probability P(t) of being cured at a given time t after diagnosis knowing
#' that he/she was alive up to time t.
#'
#' @param x result of the \code{predCuresurv} function
#'
#'
#' @param fun in "haz" or "surv" or "pt_cure", "cumhaz", "logcumhaz", the plot
#' produced is that of (excess) hazard, or that of (net) survival, or that of
#' the probability P(t) of being cured at a given time t after diagnosis
#' knowing that he/she was alive up to time t is provided, or that of
#' cumulative hazard or that of the logarithm of the cumulative hazard; if
#'  \code{fun = "all"}, the plots of the three first indicators are produced.
#'
#' @param conf.int an argument expected to be TRUE if the confidence intervals of the
#' related-indicator specified by the argument "fun" are needed. The default
#' option is FALSE. Confidence intervals are not available for \code{fun="cumhaz"} and \code{fun="logcumhaz"}
#'
#' @param conf.type One of "plain", "log", "log-log".
#' The first option causes the standard intervals curve +- k *se(curve),
#' where k is determined from conf.int. The log option calculates intervals
#'  based on log(curve). The log-log option bases
#'  the intervals on the log(-log(curve)).
#'
#' @param legend.out an argument deciding the place of the legend if \code{fun="all"}.
#' The default value is TRUE and forces most of the legend on the empty bottom-right plot slot.
#' If value is FALSE, the legend will be printed entirely in each subplot.
#'
#' @param xlab label for the x-axis of the plot.
#'
#' @param ylab optional label for the y-axis of the plot. Depending to the curve
#' of interest (hazard, survival, probability of being cured at a given time t,
#'  or all),the argument must be named \code{ylab.haz, ylab.surv, ylab.ptcure}.
#'  If missing some default labels are provided depending on the curve of
#'  interest. This name can be found in the data.frame from the result of
#'  the \code{predict.curesurv} function.
#'
#' @param ylab.haz optional label for the y-axis of the plot of excess hazard
#'
#' @param ylab.surv optional label for the y-axis of the plot of net survival
#'
#' @param ylab.ptcure optional label for the y-axis of the plot of the probability
#'  \code{P(t)} of being cured at a given time t after diagnosis knowing that
#'   he/she was alive up to time t
#'
#' @param ylab.cumhaz optional label for the y-axis of the plot of cumulative excess hazard
#'
#' @param ylab.logcumhaz optional label for the y-axis of the plot of logarithm of cumulative excess hazard
#'
#' @param col.haz optional argument to specify the color of curve of the excess hazard
#'
#' @param col.surv optional argument to specify the color of curve of the net survival
#'
#' @param col.ptcure optional argument to specify the color of curve of probability
#'  \code{P(t)} of being cured at a given time t after diagnosis knowing that
#'   he/she was alive up to time \code{t}.
#'
#' @param col.cumhaz optional argument to specify the color of curve of cumulative excess hazard
#'
#' @param col.logcumhaz optional argument to specify the color of curve of the logarithm of cumulative excess hazard
#'
#' @param col.tau optional argument to specify the color of curve of time-to-null excess hazard
#'
#' @param col.ttc optional argument to specify the color of curve of time-to-cure
#'
#' @param col.p95 optional argument to specify the color for the line highlighting \eqn{\epsilon} when \eqn{P(t) \ge 1-\epsilon}
#'
#' @param col.pi optional argument to specify the color of cure proportion
#'
#' @param lty.surv stands for line types for net survival
#'
#' @param lty.haz stands for line types for excess hazard
#'
#' @param lty.ptcure stands for line types for probability P(t) of being cured
#' at a given time t after diagnosis knowing that he/she was alive up to time t.
#'
#' @param lty.cumhaz stands for line types for cumulative excess hazard
#'
#' @param lty.logcumhaz stands for line types for logarithm cumulative excess hazard
#'
#' @param lty.ic stands for line types for confidence intervals
#'
#' @param lty.pi stands for line types for cure proportion
#'
#' @param lty.tau stands for line types for time-to-null excess hazard
#'
#' @param lty.ttc stands for line types for time-to-cure
#'
#' @param lty.p95 stands for line types for the line highlighting \eqn{\epsilon} when \eqn{P(t) \ge 1-\epsilon}
#'
#' @param lwd.main line width for the main line (haz, surv, pt_cure, cumhaz, logcumhaz)
#'
#' @param lwd.sub line width for the additionnal lines (ttc, p95, tau...)
#'
#' @param lwd.ic line width for the confidence intervals lines
#'
#' @param ... additional options as in the classical plot method.
#'
#' @keywords plot.curesurv
#'
#' @return No value is returned.
#'
#' @author Juste Goungounga, Judith Breaud, Eugenie Blandin, Olayide Boussari, Valerie Jooste
#'
#' @examples
#'
#' \donttest{
#'
#' library("curesurv")
#' library("survival")
#'
#'  testiscancer$age_crmin <- (testiscancer$age- min(testiscancer$age)) /
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
#'
#' #'  #mean of age
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
#'
#'    # predictions at time 2 years and  of age
#'
#'    newdata3 <- with(testiscancer,
#'       expand.grid(event = 0,
#'       age_crmin = seq(min(testiscancer$age_crmin),max(testiscancer$age_crmin), 0.1),
#'       time_obs  = 2))
#'
#'    pred_age_val <- predict(object = fit_m1_ad_tneh, newdata = newdata3)
#'
#'  #plot of 3 indicators for mean age
#'
#'  plot(pred_agemean, fun="all")
#'
#'
#'  #plot of net survival for mean and maximum age (comparison)
#'
#' oldpar <- par(no.readonly = TRUE)
#'
#' par(mfrow = c(2, 2),
#'     cex = 1.0)
#' plot(pred_agemax$time,
#'     pred_agemax$ex_haz,
#'     type = "l",
#'     lty = 1,
#'     lwd = 2,
#'     xlab = "Time since diagnosis",
#'     ylab = "excess hazard")
#' lines(pred_agemean$time,
#'      pred_agemean$ex_haz,
#'      type = "l",
#'      lty = 2,
#'      lwd = 2)
#'
#' legend("topright",
#'       horiz = FALSE,
#'       legend = c("hE(t) age.max = 79.9", "hE(t) age.mean = 50.8"),
#'       col = c("black", "black"),
#'       lty = c(1, 2, 1, 1, 2, 2))
#' grid()
#'
#' plot(pred_agemax$time,
#'     pred_agemax$netsurv,
#'     type = "l",
#'     lty = 1,
#'     lwd = 2,
#'     ylim = c(0, 1),
#'     xlab = "Time since diagnosis",
#'     ylab = "net survival")
#' lines(pred_agemean$time,
#'      pred_agemean$netsurv,
#'      type = "l",
#'      lty = 2,
#'      lwd = 2)
#' legend("bottomleft",
#'        horiz = FALSE,
#'        legend = c("Sn(t) age.max = 79.9", "Sn(t) age.mean = 50.8"),
#'        col = c("black", "black"),
#'       lty = c(1, 2, 1, 1, 2, 2))
#' grid()
#'
#' plot(pred_agemax$time,
#'     pred_agemax$pt_cure,
#'     type = "l",
#'     lty = 1,
#'     lwd = 2,
#'     ylim = c(0, 1), xlim = c(0,30),
#'     xlab = "Time since diagnosis",
#'     ylab = "probability of being cured P(t)")
#'
#' lines(pred_agemean$time,
#'      pred_agemean$pt_cure,
#'      type = "l",
#'      lty = 2,
#'      lwd = 2)
#'
#'
#' abline(v = pred_agemean$tau[1],
#'       lty = 2,
#'       lwd = 2,
#'       col = "blue")
#' abline(v = pred_agemean$TTC[1],
#'        lty = 2,
#'        lwd = 2,
#'        col = "red")
#' abline(v = pred_agemax$tau[1],
#'        lty = 1,
#'        lwd = 2,
#'        col = "blue")
#' abline(v = pred_agemax$TTC[1],
#'        lty = 1,
#'        lwd = 2,
#'       col = "red")
#' grid()
#'
#' legend("bottomright",
#'        horiz = FALSE,
#'        legend = c("P(t) age.max = 79.9",
#'                  "P(t) age.mean = 50.8",
#'                  "TNEH age.max = 79.9",
#'                  "TTC age.max = 79.9",
#'                  "TNEH age.mean = 50.8",
#'                  "TTC age.mean = 50.8"),
#'       col = c("black", "black", "blue", "red", "blue", "red"),
#'       lty = c(1, 2, 1, 1, 2, 2))
#'
#'
#'  val_age <- seq(min(testiscancer$age_crmin),
#'                 max(testiscancer$age_crmin), 0.1) * sd(testiscancer$age) +
#'                 min(testiscancer$age)
#'
#'
#'  pred_age_val <- predict(object = fit_m1_ad_tneh, newdata = newdata3)
#'
#'
#' par(mfrow=c(2,2))
#'  plot(val_age,
#'      pred_age_val$ex_haz, type = "l",
#'      lty=1, lwd=2,
#'      xlab = "age",
#'      ylab = "excess hazard")
#' grid()
#'
#'  plot(val_age,
#'      pred_age_val$netsurv, type = "l", lty=1,
#'      lwd=2, xlab = "age", ylab = "net survival")
#'      grid()
#'
#'  plot(val_age,
#'      pred_age_val$pt_cure, type = "l", lty=1, lwd=2,
#'      xlab = "age",
#'      ylab = "P(t)")
#'      grid()
#' par(oldpar)
#'  }
#'
#' @seealso [curesurv::predict.curesurv()], [curesurv::print.curesurv()], [curesurv::curesurv()], `browseVignettes("curesurv")`
#'
#' @import graphics
#' @export


plot.predCuresurv <- function(x,
                          fun = "all",
                          conf.int = FALSE,
                          conf.type = c("log", "log-log", "plain"),
                          legend.out=TRUE,
                          xlab = "Time since diagnosis",
                          ylab.haz = "excess hazard",
                          ylab.surv = "net survival",
                          ylab.ptcure = "P(t)",
                          ylab.cumhaz="cumulative excess hazard",
                          ylab.logcumhaz="logarithm of cumulative excess hazard",
                          col.haz = "black",
                          col.surv = "black",
                          col.ptcure = "black",
                          col.cumhaz="black",
                          col.logcumhaz="black",
                          col.tau = "red",
                          col.ttc = "green4",
                          col.p95 = "black",
                          col.pi = "blue",
                          lty.surv = 1, lty.haz = 1, lty.ptcure = 1, lty.cumhaz=1,lty.logcumhaz=1,
                          lty.pi = 2,
                          lty.tau = 2, lty.ttc = 3, lty.p95 = 4, lty.ic=5,
                          lwd.main=1,lwd.sub=1,lwd.ic=1,
                          ...) {
  if (!inherits(x, "predCuresurv"))
    stop("Primary argument much be a predCuresurv object")

  my_args <- list(...)

  if ((fun == "all") | is.null(fun)) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    if (conf.int == FALSE) {
      par(mfrow = c(2, 2))
      plot(x$time, x$ex_haz, type = "l",
           xlab = xlab,
           ylab = ylab.haz,
           col = col.haz,
           lty = lty.haz,
           lwd=lwd.main,...)
      abline(v = x$time_to_cure_ttc[1], col = col.ttc, lty = lty.ttc,lwd=lwd.sub)
      if(legend.out){
        if (x$dist[1] == "tneh") {
          abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)
        }
        legend("topright", legend = c("hazard"),col = c(col.haz),lty = c(lty.haz),lwd=lwd.main)
      }else{
        if (x$dist[1] == "tneh") {
          abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)
          legend("topright", legend = c("hazard", paste0("time=TNEH=",round(x$tau[1],2)), paste0("time=TTC=", round(x$time_to_cure_ttc[1],2))),col = c(col.haz, col.tau, col.ttc),lty = c(lty.haz, lty.tau, lty.ttc),lwd=c(lwd.main,lwd.sub,lwd.sub))
        } else if (x$model[1] != "tneh") {
          legend("topright", legend = c("hazard",  paste0("time=TTC=", round(x$time_to_cure_ttc[1],2))),col = c(col.haz,col.ttc),lty = c(lty.haz, lty.ttc),lwd=c(lwd.main,lwd.sub))
        }
      }


      plot(x$time, x$netsurv, type = "l",
           xlab = xlab,
           ylab = ylab.surv,
           col = col.surv,
           ylim = c(0,1),
           lty = lty.surv,
           lwd=lwd.main,
           ...)
      abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)
      abline(v = x$time_to_cure_ttc[1], col = col.ttc, lty = lty.ttc,lwd=lwd.sub)
      if(legend.out){
        if (x$dist[1] == "tneh") {
          abline(h = x$netsurv_tau[1], col = col.pi, lty = lty.pi,lwd=lwd.sub)
        }else if (x$model[1] != "tneh") {
          abline(h=x$cured,col=col.pi,lty=lty.pi,lwd=lwd.sub)
          }
        legend("bottomleft", legend = c("net survival"),col = c(col.surv),lty = c(lty.surv),lwd=lwd.main)

      }else{
        if (x$dist[1] == "tneh") {
          abline(h = x$netsurv_tau[1], col = col.pi, lty = lty.pi,lwd=lwd.sub)
          legend("bottomleft", legend = c("net survival",paste0("time=TNEH=",round(x$tau[1],2)),paste0("time=TTC=", round(x$time_to_cure_ttc[1],2)),paste0("cure proportion=",round(x$cured[1]*100,2),"%")),col = c(col.surv, col.tau, col.ttc, col.pi),lty = c(lty.surv, lty.tau, lty.ttc, lty.pi),lwd=c(lwd.main,rep(lwd.sub,3)))
        } else if (x$model[1] != "tneh") {
          abline(h=x$cured,col=col.pi,lty=lty.pi,lwd=lwd.sub)
          legend("bottomleft", legend = c("net survival",paste0("time=TTC=", round(x$time_to_cure_ttc[1],2)),paste0("cure proportion=",round(x$cured[1]*100,2),"%")),col = c(col.surv, col.ttc, col.pi),lty = c(lty.surv, lty.ttc, lty.pi),lwd=c(lwd.main,rep(lwd.sub,2)))
        }
      }




      plot(x$time,
           x$pt_cure, type = "l",
           xlab = xlab,
           ylab = ylab.ptcure,
           ylim = c(0, max(x$pt_cure)),
           xlim = c(0, max(x$time)),
           col = col.ptcure,
           lty = lty.ptcure, lwd=lwd.main,
           ...)
      abline(h = c(1 - x$epsilon[1]), col = col.p95, lty = lty.p95,lwd=lwd.sub)
      abline(v = x$time_to_cure_ttc[1], col = col.ttc, lty = lty.ttc,lwd=lwd.sub)
      if(legend.out){
        if (x$dist[1] == "tneh") {
          abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)}

          legend("bottomright", legend = c("P(t)"),col = c(col.ptcure),lty = c(lty.ptcure),lwd=lwd.main)

      }else{
        if (x$dist[1] == "tneh") {
          abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)
          legend("bottomright",legend = c("P(t)", paste0("time=TNEH=",round(x$tau[1],2)), paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)), paste0("P(t)=", (1 - x$epsilon[1]))),col = c(col.ptcure, col.tau, col.ttc, col.p95),lty = c(lty.ptcure, lty.tau, lty.ttc, lty.p95),lwd=c(lwd=lwd.main,rep(lwd=lwd.sub,3)))
        }else {
          legend("bottomright",legend = c("P(t)", paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)), paste0("P(t)=", (1 - x$epsilon[1]))),col = c(col.ptcure, col.ttc, col.p95),lty = c(lty.ptcure, lty.ttc, lty.p95),lwd=c(lwd=lwd.main,rep(lwd=lwd.sub,2)))

        }

      }

      if(legend.out){
      plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE,ann=FALSE)
      if (x$dist[1] == "tneh") {
        legend("center",legend = c(paste0("time=TNEH=",round(x$tau[1],2)), paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)),paste0("cure proportion=",round(x$cured[1]*100,2),"%"), paste0("P(t)=", (1 - x$epsilon[1]))),col = c( col.tau, col.ttc,col.pi, col.p95),lty = c( lty.tau, lty.ttc,lty.pi, lty.p95),lwd=rep(lwd.sub,4))
      }else {
        legend("center",legend = c( paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)),paste0("cure proportion=",round(x$cured[1]*100,2),"%"), paste0("P(t)=", (1 - x$epsilon[1]))),col = c( col.ttc,col.pi, col.p95),lty = c(lty.ttc,lty.pi, lty.p95),lwd=rep(lwd.sub,3))

      }
      }
    }
    else {

      par(mfrow = c(2, 2))
      if (x$dist[1] == "tneh") {
      plot(x$time, x$ex_haz, type = "l",
           xlab = xlab,
           ylab = ylab.haz,
           col = col.haz,
           lty = lty.haz,
           lwd=lwd.main,
           ylim=range(x$lower_bound_exhaz_tneh,x$upper_bound_exhaz_tneh),...)
      }else{
        plot(x$time, x$ex_haz, type = "l",
             xlab = xlab,
             ylab = ylab.haz,
             col = col.haz,
             lty = lty.haz,
             lwd=lwd.main,
             ylim=range(x$lower_bound_exhaz_ic,x$upper_bound_exhaz_ic),...)
      }
      abline(v = x$time_to_cure_ttc[1], col = col.ttc, lty = lty.ttc,lwd=lwd.sub)
      if(legend.out){
        if (x$dist[1] == "tneh") {
          abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)
          if(conf.type=="plain"){
            lines(x$time,x$lower_bound_exhaz_tneh,col=col.haz,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_exhaz_tneh,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log"){
            lines(x$time,x$lower_bound_exhaz_tneh_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_exhaz_tneh_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log-log"){
            lines(x$time,x$lower_bound_exhaz_tneh_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_exhaz_tneh_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          }
        }else if (x$model[1] != "tneh"){
          if(conf.type=="plain"){
            lines(x$time,x$lower_bound_exhaz_ic,col=col.haz,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_exhaz_ic,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log"){
            lines(x$time,x$lower_bound_exhaz_ic_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_exhaz_ic_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log-log"){
            lines(x$time,x$lower_bound_exhaz_ic_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_exhaz_ic_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          }
        }
        legend("topright", legend = c("hazard","IC"),col = c(col.haz,col.haz),lty = c(lty.haz,lty.ic),lwd=c(lwd.main,lwd.ic))
      }else{
        if (x$dist[1] == "tneh") {
          abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)
          if(conf.type=="plain"){
            lines(x$time,x$lower_bound_exhaz_tneh,col=col.haz,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_exhaz_tneh,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log"){
            lines(x$time,x$lower_bound_exhaz_tneh_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_exhaz_tneh_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log-log"){
            lines(x$time,x$lower_bound_exhaz_tneh_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_exhaz_tneh_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          }
          legend("topright", legend = c("hazard","IC", paste0("time=TNEH=",round(x$tau[1],2)), paste0("time=TTC=",round(x$time_to_cure_ttc[1],2))),col = c(col.haz,col.haz, col.tau, col.ttc),lty = c(lty.haz,lty.ic, lty.tau, lty.ttc),lwd=c(lwd.main,lwd.ic,lwd.sub,lwd.sub))
        } else if (x$model[1] != "tneh") {
          if(conf.type=="plain"){
            lines(x$time,x$lower_bound_exhaz_ic,col=col.haz,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_exhaz_ic,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log"){
            lines(x$time,x$lower_bound_exhaz_ic_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_exhaz_ic_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log-log"){
            lines(x$time,x$lower_bound_exhaz_ic_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_exhaz_ic_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          }
          legend("topright", legend = c("hazard","IC",  paste0("time=TTC=",round(x$time_to_cure_ttc[1],2))),col = c(col.haz,col.haz,col.ttc),lty = c(lty.haz,lty.ic, lty.ttc),lwd=c(lwd.main,lwd.ic,lwd.sub))
        }
      }


      plot(x$time, x$netsurv, type = "l",
           xlab = xlab,
           ylab = ylab.surv,
           col = col.surv,
           ylim = c(0,1),
           lty = lty.surv,
           lwd=lwd.main,
           ...)
      abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)
      abline(v = x$time_to_cure_ttc[1], col = col.ttc, lty = lty.ttc,lwd=lwd.sub)
      if(legend.out){
        if (x$dist[1] == "tneh") {
          if(conf.type=="plain"){
            lines(x$time,x$lower_bound_sn_tneh,col=col.surv,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_sn_tneh,col=col.surv,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log"){
            lines(x$time,x$lower_bound_sn_tneh_log,col=col.surv,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_sn_tneh_log,col=col.surv,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log-log"){
            lines(x$time,x$lower_bound_sn_tneh_loglog,col=col.surv,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_sn_tneh_loglog,col=col.surv,lty=lty.ic,lwd=lwd.ic)
          }
          abline(h = x$netsurv_tau[1], col = col.pi, lty = lty.pi,lwd=lwd.sub)
        }else if (x$model[1] != "tneh") {
          abline(h=x$cured,col=col.pi,lty=lty.pi,lwd=lwd.sub)
          if(conf.type=="plain"){
            lines(x$time,x$lower_bound_sn_ic,col=col.surv,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_sn_ic,col=col.surv,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log"){
            lines(x$time,x$lower_bound_sn_ic_log,col=col.surv,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_sn_ic_log,col=col.surv,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log-log"){
            lines(x$time,x$lower_bound_sn_ic_loglog,col=col.surv,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_sn_ic_loglog,col=col.surv,lty=lty.ic,lwd=lwd.ic)
          }
        }
        legend("bottomleft", legend = c("net survival","IC"),col = rep(col.surv,2),lty = c(lty.surv,lty.ic),lwd=c(lwd.main,lwd.ic))

      }else{
        if (x$dist[1] == "tneh") {
          abline(h = x$netsurv_tau[1], col = col.pi, lty = lty.pi,lwd=lwd.sub)
          if(conf.type=="plain"){
            lines(x$time,x$lower_bound_sn_tneh,col=col.surv,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_sn_tneh,col=col.surv,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log"){
            lines(x$time,x$lower_bound_sn_tneh_log,col=col.surv,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_sn_tneh_log,col=col.surv,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log-log"){
            lines(x$time,x$lower_bound_sn_tneh_loglog,col=col.surv,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_sn_tneh_loglog,col=col.surv,lty=lty.ic,lwd=lwd.ic)
          }
          legend("bottomleft", legend = c("net survival","IC",paste0("time=TNEH=",round(x$tau[1],2)),paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)),paste0("cure proportion=",round(x$cured[1]*100,2),"%")),col = c(col.surv,col.surv, col.tau, col.ttc, col.pi),lty = c(lty.surv,lty.ic, lty.tau, lty.ttc, lty.pi),lwd=c(lwd.main,lwd.ic,rep(lwd.sub,3)))
        } else if (x$model[1] != "tneh") {
          abline(h=x$cured,col=col.pi,lty=lty.pi,lwd=lwd.sub)
          if(conf.type=="plain"){
            lines(x$time,x$lower_bound_sn_ic,col=col.surv,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_sn_ic,col=col.surv,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log"){
            lines(x$time,x$lower_bound_sn_ic_log,col=col.surv,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_sn_ic_log,col=col.surv,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log-log"){
            lines(x$time,x$lower_bound_sn_ic_loglog,col=col.surv,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_sn_ic_loglog,col=col.surv,lty=lty.ic,lwd=lwd.ic)
          }
          legend("bottomleft", legend = c("net survival","IC",paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)),paste0("cure proportion=",round(x$cured[1]*100,2),"%")),col = c(col.surv,col.surv, col.ttc, col.pi),lty = c(lty.surv,lty.ic, lty.ttc, lty.pi),lwd=c(lwd.main,lwd.ic,rep(lwd.sub,2)))
        }
      }




      plot(x$time,
           x$pt_cure, type = "l",
           xlab = xlab,
           ylab = ylab.ptcure,
           ylim = c(0, max(x$pt_cure)),
           xlim = c(0, max(x$time)),
           col = col.ptcure,
           lty = lty.ptcure, lwd=lwd.main,
           ...)
      abline(h = c(1 - x$epsilon[1]), col = col.p95, lty = lty.p95,lwd=lwd.sub)
      abline(v = x$time_to_cure_ttc[1], col = col.ttc, lty = lty.ttc,lwd=lwd.sub)
      if(legend.out){
        if (x$dist[1] == "tneh") {
          abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)
        if(conf.type=="plain"){
          lines(x$time,x$lower_bound_pt_cure_tneh,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_pt_cure_tneh,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
        }else if(conf.type=="log"){
          lines(x$time,x$lower_bound_pt_cure_tneh_log,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_pt_cure_tneh_log,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
        }else if(conf.type=="log-log"){
          lines(x$time,x$lower_bound_pt_cure_tneh_loglog,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_pt_cure_tneh_loglog,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
        }}else if (x$model[1] != "tneh") {
          if(conf.type=="plain"){
            lines(x$time,x$lower_bound_pt_ic,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_pt_ic,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log"){
            lines(x$time,x$lower_bound_pt_ic_log,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_pt_ic_log,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log-log"){
            lines(x$time,x$lower_bound_pt_ic_loglog,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_pt_ic_loglog,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          }
        }

        legend("bottomright", legend = c("P(t)","IC"),col = rep(col.ptcure,2),lty = c(lty.ptcure,lty.ic),lwd=c(lwd.main,lwd.ic))

      }else{
        if (x$dist[1] == "tneh") {
          abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)
          if(conf.type=="plain"){
            lines(x$time,x$lower_bound_pt_cure_tneh,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_pt_cure_tneh,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log"){
            lines(x$time,x$lower_bound_pt_cure_tneh_log,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_pt_cure_tneh_log,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log-log"){
            lines(x$time,x$lower_bound_pt_cure_tneh_loglog,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_pt_cure_tneh_loglog,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          }
          legend("bottomright",legend = c("P(t)","IC", paste0("time=TNEH=",round(x$tau[1],2)), paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)), paste0("P(t)=", (1 - x$epsilon[1]))),col = c(col.ptcure,col.ptcure, col.tau, col.ttc, col.p95),lty = c(lty.ptcure,lty.ic, lty.tau, lty.ttc, lty.p95),lwd=c(lwd=lwd.main,lwd.ic,rep(lwd=lwd.sub,3)))
        }else {
          if(conf.type=="plain"){
            lines(x$time,x$lower_bound_pt_ic,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_pt_ic,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log"){
            lines(x$time,x$lower_bound_pt_ic_log,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_pt_ic_log,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          }else if(conf.type=="log-log"){
            lines(x$time,x$lower_bound_pt_ic_loglog,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
            lines(x$time,x$upper_bound_pt_ic_loglog,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          }
          legend("bottomright",
                 legend = c("P(t)","IC", paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)), paste0("P(t)=", (1 - x$epsilon[1]))),
                 col = c(col.ptcure,col.ptcure, col.ttc, col.p95),
                 lty = c(lty.ptcure,lty.ic, lty.ttc, lty.p95),
                 lwd=c(lwd=lwd.main,lwd.ic,rep(lwd=lwd.sub,2)))

        }

      }

      if(legend.out){
        plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE,ann=FALSE)
        if (x$dist[1] == "tneh") {
          legend("center",legend = c(paste0("time=TNEH=",round(x$tau[1],2)), paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)),paste0("cure proportion=",round(x$cured[1]*100,2),"%"), paste0("P(t)=", (1 - x$epsilon[1]))),col = c( col.tau, col.ttc,col.pi, col.p95),lty = c( lty.tau, lty.ttc,lty.pi, lty.p95),lwd=rep(lwd.sub,4))
        }else {
          legend("center",legend = c( paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)),paste0("cure proportion=",round(x$cured[1]*100,2),"%"), paste0("P(t)=", (1 - x$epsilon[1]))),col = c( col.ttc,col.pi, col.p95),lty = c(lty.ttc,lty.pi, lty.p95),lwd=rep(lwd.sub,3))

        }
      }
    }


  } else if (fun == "haz")
    {

    if (conf.int == FALSE) {
      plot(x$time,
           x$ex_haz, type = "l",
           xlab = xlab, ylab = ylab.haz,
           col = col.haz, lty = lty.haz, lwd=lwd.main, ...)
      abline(v = x$time_to_cure_ttc[1], col = col.ttc, lty = 2,lwd=lwd.sub)

      if (x$dist[1] == "tneh") {
        abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)
        legend("topright", legend = c("hazard", paste0("time=TNEH=",round(x$tau[1],2)),
                                      paste0("time=TTC=",round(x$time_to_cure_ttc[1],2))),
               col = c(col.haz, col.tau, col.ttc),
               lty = c(lty.haz, lty.tau,lty.ttc),
               lwd=c(lwd.main,rep(lwd.sub,2)))
      } else if (x$model[1] != "tneh") {
        legend("topright", legend = c("hazard",  paste0("time=TTC=",round(x$time_to_cure_ttc[1],2))),
               col = c(col.haz,col.ttc),
               lty = c(lty.haz, lty.ttc),
               lwd=c(lwd.main,lwd.sub))
      }


    } else {


      if (x$dist[1] == "tneh") {
        plot(x$time,
             x$ex_haz, type = "l",
             xlab = xlab, ylab = ylab.haz,
             col = col.haz, lty = lty.haz, lwd=lwd.main,
             ylim=range(x$lower_bound_exhaz_tneh,x$upper_bound_exhaz_tneh), ...)
        abline(v = x$time_to_cure_ttc[1], col = col.ttc, lty = 2,lwd=lwd.sub)
        abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)
        if(conf.type=="plain"){
          lines(x$time,x$lower_bound_exhaz_tneh,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_exhaz_tneh,col=col.haz,lty=lty.ic,lwd=lwd.ic)
        }else if(conf.type=="log"){
          lines(x$time,x$lower_bound_exhaz_tneh_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_exhaz_tneh_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
        }else if(conf.type=="log-log"){
          lines(x$time,x$lower_bound_exhaz_tneh_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_exhaz_tneh_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
        }
        legend("topright", legend = c("hazard","IC hazard", paste0("time=TNEH=",round(x$tau[1],2)),
                                      paste0("time=TTC=",round(x$time_to_cure_ttc[1],2))),
               col = c(col.haz,col.haz, col.tau, col.ttc),
               lty = c(lty.haz,lty.ic, lty.tau,lty.ttc),
               lwd=c(lwd.main,lwd.ic,rep(lwd.sub,2)))
      } else if (x$model[1] != "tneh") {
        plot(x$time,
             x$ex_haz, type = "l",
             xlab = xlab, ylab = ylab.haz,
             col = col.haz, lty = lty.haz, lwd=lwd.main,
             ylim=range(x$lower_bound_exhaz_ic,x$upper_bound_exhaz_ic), ...)
        abline(v = x$time_to_cure_ttc[1], col = col.ttc, lty = 2,lwd=lwd.sub)

        if(conf.type=="plain"){
          lines(x$time,x$lower_bound_exhaz_ic,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_exhaz_ic,col=col.haz,lty=lty.ic,lwd=lwd.ic)
        }else if(conf.type=="log"){
          lines(x$time,x$lower_bound_exhaz_ic_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_exhaz_ic_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
        }else if(conf.type=="log-log"){
          lines(x$time,x$lower_bound_exhaz_ic_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_exhaz_ic_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
        }
        legend("topright", legend = c("hazard","IC hazard",  paste0("time=TTC=",round(x$time_to_cure_ttc[1],2))),
               col = c(col.haz,col.haz,col.ttc),
               lty = c(lty.haz,lty.ic, lty.ttc),
               lwd=c(lwd.main,lwd.ic,lwd.sub))
      }



    }

  } else if (fun == "surv")
    {
    if (conf.int == FALSE) {
      plot(x$time, x$netsurv, type = "l",
           xlab = xlab,
           ylab = ylab.surv,
           col = col.surv,
           ylim = c(0,1), lty = lty.surv,lwd=lwd.main, ...)
      abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)
      abline(v = x$time_to_cure_ttc[1], col = col.ttc, lty = lty.ttc,lwd=lwd.sub)

      if (x$dist[1] == "tneh") {
        abline(h = x$netsurv_tau[1], col = col.pi, lty = lty.pi,lwd=lwd.sub)
        legend("bottomleft", legend = c("net survival",
                                        paste0("time=TNEH=",round(x$tau[1],2)),
                                        paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)),
                                      paste0("cure proportion=",round(x$cured[1]*100,2),"%")),
               col = c(col.surv, col.tau, col.ttc, col.pi),
               lty = c(lty.surv, lty.tau, lty.ttc, lty.pi),
               lwd=c(lwd.main,rep(lwd.sub,3)))
      } else if (x$model[1] != "tneh")
        {
        abline(h=x$cured,col=col.pi,lty=lty.pi,lwd=lwd.sub)

        legend("bottomleft", legend = c("net survival",
                                        paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)),
                                      paste0("cure proportion=",round(x$cured[1]*100,2),"%")),
               col = c(col.surv, col.ttc, col.pi),
               lty = c(lty.surv, lty.ttc, lty.pi),
               lwd=c(lwd.main,rep(lwd.sub,2)))
      }


    }
    else {
      plot(x$time, x$netsurv, type = "l",
           xlab = xlab,
           ylab = ylab.surv,
           col = col.surv,
           ylim = c(0,1), lty = lty.surv, lwd=lwd.main, ...)
      abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)
      abline(v = x$time_to_cure_ttc[1], col = col.ttc, lty = lty.ttc,lwd=lwd.sub)

      if (x$dist[1] == "tneh") {
        abline(h = x$netsurv_tau[1], col = col.pi, lty = lty.pi,lwd=lwd.sub)
        if(conf.type=="plain"){
          lines(x$time,x$lower_bound_sn_tneh,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_sn_tneh,col=col.haz,lty=lty.ic,lwd=lwd.ic)
        }
        else if(conf.type=="log"){
          lines(x$time,x$lower_bound_sn_tneh_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_sn_tneh_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
        } else if(conf.type=="log-log")
          {
          lines(x$time,x$lower_bound_sn_tneh_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_sn_tneh_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
        }else{
          stop("not implemented for this conf.type")
        }

        legend("bottomleft", legend = c("net survival", "IC net suvival",
                                        paste0("time=TNEH=",round(x$tau[1],2)),
                                        paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)),
                                        paste0("cure proportion=",round(x$cured[1]*100,2),"%")),
               col = c(col.surv,col.surv, col.tau, col.ttc, col.pi),
               lty = c(lty.surv,lty.ic, lty.tau, lty.ttc, lty.pi),
               lwd=c(lwd.main,lwd.ic,rep(lwd.sub,3)))
      }
      else if (x$model[1] != "tneh") {
        abline(h=x$cured,col=col.pi,lty=lty.pi,lwd=lwd.sub)
        if(conf.type=="plain"){
          lines(x$time,x$lower_bound_sn_ic,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_sn_ic,col=col.haz,lty=lty.ic,lwd=lwd.ic)
        }
        else if(conf.type=="log"){
          lines(x$time,x$lower_bound_sn_ic_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_sn_ic_log,col=col.haz,lty=lty.ic,lwd=lwd.ic)
        } else if(conf.type=="log-log")
          {
          lines(x$time,x$lower_bound_sn_ic_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_sn_ic_loglog,col=col.haz,lty=lty.ic,lwd=lwd.ic)
        }else
          {
          stop("not implemented for this conf.type")
        }

        legend("bottomleft", legend = c("net survival", "IC net survival",
                                        paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)),
                                        paste0("cure proportion=",round(x$cured[1]*100,2),"%")),
               col = c(col.surv,col.surv, col.ttc, col.pi),
               lty = c(lty.surv,lty.ic, lty.ttc, lty.pi),
               lwd=c(lwd.main,lwd.ic,rep(lwd.sub,2)))
      }
      }
    } else if (fun == "pt_cure")
    {

    if (conf.int == FALSE) {
      plot(x$time,
           x$pt_cure, type = "l",
           xlab = xlab,
           ylab = ylab.ptcure,
           ylim = c(0, max(x$pt_cure)),
           xlim = c(0, max(x$time)),
           col = col.ptcure,
           lty = lty.ptcure, lwd=lwd.main, ...)
      abline(h = c(1 - x$epsilon[1]), col = col.p95, lty = lty.p95,lwd=lwd.sub)
      abline(v = x$time_to_cure_ttc[1], col = col.ttc, lty = lty.ttc,lwd=lwd.sub)

      if (x$dist[1] == "tneh") {
        abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)

        legend("bottomright",
               legend = c("P(t)", paste0("time=TNEH=",round(x$tau[1],2)),
                          paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)), paste0("P(t)=", (1 - x$epsilon[1]))),
               col = c(col.ptcure, col.tau, col.ttc, col.p95),
               lty = c(lty.ptcure,lty.tau, lty.ttc, lty.p95),
               lwd=c(lwd.main,rep(lwd.sub,3)))
      }
      else {

        legend("bottomright",
               legend = c("P(t)", paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)), paste0("P(t)=", (1 - x$epsilon[1]))),
               col = c(col.ptcure, col.ttc, col.p95),
               lty = c(lty.ptcure, lty.ttc, lty.p95),
               lwd=c(lwd.main,rep(lwd.sub,2)))
      }



    }
      else {
      plot(x$time,
           x$pt_cure, type = "l",
           xlab = xlab,
           ylab = ylab.ptcure,
           ylim = c(0, max(x$pt_cure)),
           xlim = c(0, max(x$time)),
           col = col.ptcure,
           lty = lty.ptcure, lwd=lwd.main, ...)
      abline(h = c(1 - x$epsilon[1]), col = col.p95, lty = lty.p95,lwd=lwd.sub)
      abline(v = x$time_to_cure_ttc[1], col = col.ttc, lty = lty.ttc,lwd=lwd.sub)

      if (x$dist[1] == "tneh") {
        abline(v = x$tau[1], col = col.tau, lty = lty.tau,lwd=lwd.sub)
        if(conf.type=="plain"){
          lines(x$time,x$lower_bound_pt_cure_tneh,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_pt_cure_tneh,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
        }
        else if(conf.type=="log"){
          lines(x$time,x$lower_bound_pt_cure_tneh_log,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_pt_cure_tneh_log,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
        } else if(conf.type=="log-log")
        {
          lines(x$time,x$lower_bound_pt_cure_tneh_loglog,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_pt_cure_tneh_loglog,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
        }else{
          stop("not implemented for this conf.type")
        }


        legend("bottomright",
               legend = c("P(t)","IC P(t)", paste0("time=TNEH=",round(x$tau[1],2)),
                          paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)), paste0("P(t)=", (1 - x$epsilon[1]))),
               col = c(col.ptcure,col.ptcure, col.tau, col.ttc, col.p95),
               lty = c(lty.ptcure,lty.ic,lty.tau, lty.ttc, lty.p95),
               lwd=c(lwd.main,lwd.ic,rep(lwd.sub,3)))
      }else if (x$model[1] != "tneh") {
        if(conf.type=="plain"){
          lines(x$time,x$lower_bound_pt_ic,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_pt_ic,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
        }
        else if(conf.type=="log"){
          lines(x$time,x$lower_bound_pt_ic_log,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_pt_ic_log,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
        } else if(conf.type=="log-log")
        {
          lines(x$time,x$lower_bound_pt_ic_loglog,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_pt_ic_loglog,col=col.ptcure,lty=lty.ic,lwd=lwd.ic)
        }else{
          stop("not implemented for this conf.type")
        }
        legend("bottomright",
               legend = c("P(t)","IC P(t)", paste0("time=TTC=",round(x$time_to_cure_ttc[1],2)), paste0("P(t)=", (1 - x$epsilon[1]))),
               col = c(col.ptcure,col.ptcure, col.ttc, col.p95),
               lty = c(lty.ptcure,lty.ic, lty.ttc, lty.p95),
               lwd=c(lwd.main,lwd.ic,rep(lwd.sub,2)))
      }
    }

  } else if (fun == "cumhaz")
    {

    if (conf.int == FALSE) {
      plot(x$time, -log(x$netsurv), type = "l",
           xlab = xlab,
           ylab = ylab.cumhaz,
           col = col.cumhaz,lty=lty.cumhaz,lwd=lwd.main, ...)
      legend("bottomright",col=col.cumhaz,lty=lty.cumhaz,legend="cumulative excess hazard",lwd=lwd.main)


    }
    else {
      stop("confidence intervals not yet implemented for fun='cumhaz'")
    }

  } else if (fun == "logcumhaz")
    {
    if (conf.int == FALSE) {
      plot(x$time, log(-log(x$netsurv)), type = "l",
           xlab = xlab,
           ylab = ylab.logcumhaz,
           col = col.logcumhaz,lty=lty.logcumhaz,lwd=lwd.main,...)
      legend("bottomright",legend=c("log of cumulative excess hazard"),col=c(col.logcumhaz),lty=c(lty.logcumhaz),lwd=lwd.main)
    }
    else {
      if(conf.int==TRUE & conf.type=="plain"){
        if (x$dist[1] == "tneh") {
          stop("confidence intervals are not yet implemented for the tneh model when fun='logcumhaz'")
        }else if (x$model[1] != "tneh"){
          plot(x$time, log(-log(x$netsurv)), type = "l",
               xlab = xlab,
               ylab = ylab.logcumhaz,
               col = col.logcumhaz,lty=lty.logcumhaz,lwd=lwd.main,...)
          lines(x$time,x$lower_bound_logCumHaz_ic,col=col.logcumhaz,lty=lty.ic,lwd=lwd.ic)
          lines(x$time,x$upper_bound_logCumHaz_ic,col=col.logcumhaz,lty=lty.ic,lwd=lwd.ic)
          legend("bottomright",legend=c("log of cumulative excess hazard","IC"),col=c(col.logcumhaz,col.logcumhaz),lty=c(lty.logcumhaz,lty.ic),lwd=c(lwd.main,lwd.ic))
        }

      }else{
        stop("only conf.type='plain' is implemented for fun='logcumhaz'")
      }

    }
  }
  invisible()

}

