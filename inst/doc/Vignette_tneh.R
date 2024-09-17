## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(curesurv)

## -----------------------------------------------------------------------------
fit_ad_tneh_nocov <- curesurv(Surv(time_obs, event) ~ 1,
                             pophaz = "ehazard",
                             cumpophaz = "cumehazard",
                             model = "nmixture", dist = "tneh",  
                             link_tau = "linear",
                             data = testiscancer,
                             method_opt = "L-BFGS-B")



## -----------------------------------------------------------------------------
summary(fit_ad_tneh_nocov)

## -----------------------------------------------------------------------------
newdata1 <- with(testiscancer,
  expand.grid(event = 0, time_obs  = seq(0.1,10,0.1)))
p_28 <- predict(object = fit_ad_tneh_nocov, newdata = newdata1)


## -----------------------------------------------------------------------------
plot(p_28)

## -----------------------------------------------------------------------------
p_28[1,c("lower_bound_pi_tneh","cured","upper_bound_pi_tneh")]

## -----------------------------------------------------------------------------
p_28[1,c("lower_bound_TTC_tneh","time_to_cure_ttc","upper_bound_TTC_tneh")]

## -----------------------------------------------------------------------------
 testiscancer$age_crmin <- (testiscancer$age- min(testiscancer$age)) /
              sd(testiscancer$age)

## -----------------------------------------------------------------------------
oldpar<-par(no.readonly = TRUE)
par(mfrow=c(1,1))
plot(testiscancer$age,testiscancer$age_crmin,type="l",xlab="age",ylab="age_crmin")
text(x=min(testiscancer$age),y=0,col="blue",round(min(testiscancer$age),2))
lines(x=rep(min(testiscancer$age),2),y=c(-1,0),col="blue",lty=2)
lines(x=c(-1,min(testiscancer$age)),y=c(0,0),col="blue",lty=2)
text(x=min(testiscancer$age)+sd(testiscancer$age),y=1.2,col="red",round(sd(testiscancer$age)+min(testiscancer$age),2))
lines(x=rep(min(testiscancer$age)+sd(testiscancer$age),2),y=c(-1,1),col="red",lty=2)
lines(x=c(-1,min(testiscancer$age)+sd(testiscancer$age)),y=c(1,1),col="red",lty=2)
par(oldpar)

## -----------------------------------------------------------------------------
fit_m1_ad_tneh <- curesurv(Surv(time_obs, event) ~ z_tau(age_crmin) + 
                          z_alpha(age_crmin),
                          pophaz = "ehazard",
                          cumpophaz = "cumehazard",
                          model = "nmixture", dist = "tneh",
                          link_tau = "linear",
                          data = testiscancer,
                          method_opt = "L-BFGS-B")

## -----------------------------------------------------------------------------
 summary(fit_m1_ad_tneh)

## -----------------------------------------------------------------------------
  #time varying prediction for patient with age mean
newdata1 <- with(testiscancer, 
                 expand.grid(event = 0,
                             age_crmin = mean(age_crmin),
                             time_obs  = seq(0.001,10,0.1)))

 pred_agemean <- predict(object = fit_m1_ad_tneh, newdata = newdata1)

## -----------------------------------------------------------------------------
#time varying prediction for patient with age max
newdata2 <- with(testiscancer, 
                  expand.grid(event = 0,
                              age_crmin = max(age_crmin),
                              time_obs  = seq(0.001,10,0.1)))

pred_agemax <- predict(object = fit_m1_ad_tneh, newdata = newdata2)

## -----------------------------------------------------------------------------
# predictions at time 2 years with varying age

   newdata3 <- with(testiscancer,
      expand.grid(event = 0, 
                  age_crmin = seq(min(testiscancer$age_crmin), 
                                  max(testiscancer$age_crmin), 0.1),
                  time_obs  = 2))

pred_age_val <- predict(object = fit_m1_ad_tneh, newdata = newdata3)
val_age <- seq(min(testiscancer$age_crmin),
               max(testiscancer$age_crmin),
               0.1) * sd(testiscancer$age) +  min(testiscancer$age)

## -----------------------------------------------------------------------------
par(mfrow = c(2, 2),cex = 1.0)

plot(pred_agemax$time,pred_agemax$ex_haz,type = "l",lty = 1,lwd = 2,
     xlab = "Time since diagnosis", ylab = "excess hazard")

lines(pred_agemean$time,pred_agemean$ex_haz,type = "l",lty = 2,lwd = 2)

legend("topright",horiz = FALSE,
       legend = c("hE(t) age.max = 79.9", "hE(t) age.mean = 50.8"),
       col = c("black", "black"),lty = c(1, 2, 1, 1, 2, 2))
grid()

plot(pred_agemax$time,pred_agemax$netsurv,type = "l",lty = 1,lwd = 2,
    ylim = c(0, 1),xlab = "Time since diagnosis",ylab = "net survival")

lines(pred_agemean$time,pred_agemean$netsurv,type = "l",lty = 2,lwd = 2)

legend("bottomleft",horiz = FALSE,
       legend = c("Sn(t) age.max = 79.9", "Sn(t) age.mean = 50.8"),
       col = c("black", "black"),lty = c(1, 2, 1, 1, 2, 2))
grid()

plot(pred_agemax$time,pred_agemax$pt_cure,type = "l",lty = 1,lwd = 2,ylim = c(0, 1),
     xlab = "Time since diagnosis",ylab = "probability of being cured P(t)")

lines(pred_agemean$time,pred_agemean$pt_cure,type = "l",lty = 2,lwd = 2)

abline(v = pred_agemean$tau[1],lty = 2,lwd = 2,col = "blue")
abline(v = pred_agemean$time_to_cure_ttc[1],lty = 2,lwd = 2,col = "red")
abline(v = pred_agemax$tau[1],lty = 1,lwd = 2,col = "blue")
abline(v = pred_agemax$time_to_cure_ttc[1],lty = 1,lwd = 2,col = "red")
grid()

legend("bottomright",horiz = FALSE,
       legend = c("P(t) age.max = 79.9","P(t) age.mean = 50.8",
                  "TNEH age.max = 79.9","TTC age.max = 79.9",
                 "TNEH age.mean = 50.8","TTC age.mean = 50.8"),
      col = c("black", "black", "blue", "red", "blue", "red"),
      lty = c(1, 2, 1, 1, 2, 2))
par(oldpar)

## -----------------------------------------------------------------------------
par(mfrow=c(2,2))
plot(val_age,pred_age_val$ex_haz, 
     type = "l",lty=1, lwd=2,
     xlab = "age",ylab = "excess hazard 2y after diagnosis")
grid()

 plot(val_age,pred_age_val$netsurv, 
      type = "l", lty=1,lwd=2,ylim=c(0,1), 
      xlab = "age", ylab = "net survival 2y after diagnosis")
     grid()

 plot(val_age,pred_age_val$pt_cure, type = "l", lty=1, lwd=2,ylim=c(0,1),
     xlab = "age",ylab = "P(t) 2y after diagnosis")
     grid()
     
plot(val_age,pred_age_val$cured, type = "l", lty=1, lwd=2,ylim=c(0,1),
     xlab = "age", ylab = "cure proportion")
     grid()
     
par(oldpar)

## -----------------------------------------------------------------------------

#| echo: true
#| label: withtauonly
#| warning: false
#| message: false

fit_ad_tneh_covtau <- curesurv(
  Surv(time_obs, event) ~ z_tau(age_cr),
  pophaz = "ehazard",
  cumpophaz = "cumehazard",
  model = "nmixture",
  dist = "tneh",
  link_tau = "linear",
  data = testiscancer,
  method_opt = "L-BFGS-B"
)



## -----------------------------------------------------------------------------
summary(fit_ad_tneh_covtau)


## -----------------------------------------------------------------------------
summary(testiscancer$age_cr)
summary(testiscancer$age)
newdata2 <- with(testiscancer,
                 expand.grid(event = 0, 
                             time_obs  = seq(0.1, 10, 0.1),
                             age_cr = c(-0.9358, 0.0276, 0.9525) ))
newdata2_1stqu <- newdata2[newdata2$age_cr==-0.9358,]
newdata2_2rdqu <- newdata2[newdata2$age_cr==0.0276,]
newdata2_3rdqu <- newdata2[newdata2$age_cr==0.9525,]

p1stqu <- predict(object = fit_ad_tneh_covtau, newdata = newdata2_1stqu)
p2rdqu <- predict(object = fit_ad_tneh_covtau, newdata = newdata2_2rdqu)
p3rdqu <- predict(object = fit_ad_tneh_covtau, newdata = newdata2_3rdqu)


## -----------------------------------------------------------------------------
par(mfrow=c(2,2))

  plot(p1stqu, 
     main = "Excess hazard for age 33.3", 
     fun = "haz")

  plot(p2rdqu,
     fun = "haz",
     main = "Excess hazard for age 51.57")

  plot(p3rdqu,
     fun = "haz",
     main = "Excess hazard for age 69.11")

par(oldpar)

## -----------------------------------------------------------------------------

#| echo: true
#| label: only_covariate_on_alpha
#| message: false
#| warning: false

fit_ad_tneh_covalpha <-
  curesurv(
    Surv(time_obs, event) ~ z_alpha(age_cr),
    pophaz = "ehazard",
    cumpophaz = "cumehazard",
    model = "nmixture",
    dist = "tneh",
    link_tau = "linear",
    data = testiscancer,
    method_opt = "L-BFGS-B"
  )



## -----------------------------------------------------------------------------
summary(fit_ad_tneh_covalpha)

## -----------------------------------------------------------------------------
p4_33.3 <- predict(object = fit_ad_tneh_covalpha,
                 newdata = newdata2_1stqu)
p4_51.6 <- predict(object = fit_ad_tneh_covalpha,
                 newdata = newdata2_2rdqu)
p4_69.1 <- predict(object = fit_ad_tneh_covalpha,
                 newdata = newdata2_3rdqu)


## -----------------------------------------------------------------------------
par(mfrow=c(2,2))

  plot(p4_33.3, 
     main = "Pt_cure for age 33.3", 
     fun = "pt_cure")

  plot(p4_51.6,
     fun = "pt_cure",
     main = "Pt_cure for age 51.6")

  plot(p4_69.1,
     fun = "pt_cure",
     main = "Pt_cure for age 69.1")
par(oldpar)

## -----------------------------------------------------------------------------
AIC(fit_ad_tneh_nocov,fit_ad_tneh_covalpha,fit_ad_tneh_covtau,fit_m1_ad_tneh)

## -----------------------------------------------------------------------------
date()
sessionInfo()

