## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("curesurv_0.1.0.tar.gz",
#                   repos = NULL,
#                   type = "source")

## -----------------------------------------------------------------------------
library(xhaz)
library(survexp.fr)
library(curesurv)

## -----------------------------------------------------------------------------
# We had these information in the dataset
data("testiscancer", package = "curesurv")

head(testiscancer)

dim(testiscancer)
testiscancer$sex <- "male"
levels(testiscancer$sex) <- c("male", "female")
testiscancer$year <- as.Date("2000-01-01")

## -----------------------------------------------------------------------------

attributes(survexp.fr)


fit.haz <- exphaz(
    formula = Surv(time_obs, event) ~ 1,
    data = testiscancer,
    ratetable = survexp.fr,
    only_ehazard = FALSE,
    rmap = list(age = 'age',
                sex = 'sex', 
                year = 'year'))

## -----------------------------------------------------------------------------
#instantaneous population hazard
testiscancer$haz <- testiscancer$ehazard
#cumulative population hazard
testiscancer$cumhaz <- testiscancer$cumehazard


## -----------------------------------------------------------------------------
fit_mod0 <- curesurv(Surv(time_obs, event) ~ age_cr | age_cr,
                   pophaz = "haz",
                   cumpophaz = "cumhaz",
                   model = "mixture", dist = "weib",
                   data = testiscancer,
                   method_opt = "L-BFGS-B")

fit_mod0


## -----------------------------------------------------------------------------
fit_mod1 <- curesurv(Surv(time_obs, event) ~ age_cr | age_cr,
                   pophaz = "haz",
                   cumpophaz = "cumhaz",
                   model = "mixture", dist = "weib",
                   data = testiscancer, 
                   pophaz.alpha = TRUE,
                   method_opt = "L-BFGS-B")

fit_mod1


## -----------------------------------------------------------------------------
AIC(fit_mod0,fit_mod1)

## -----------------------------------------------------------------------------
val_age <- seq(-1.39, 1.8, 0.1)
newage <- round(val_age * sd(testiscancer$age) + mean(testiscancer$age), 2)

newdata3 <- with(testiscancer,
                 expand.grid(
                   event = 0,
                   age_cr = val_age,
    time_obs  = 2))


pred_age_mod0 <- predict(object = fit_mod0,
                                  newdata = newdata3)


pred_age_mod1 <- predict(object = fit_mod1,
                                  newdata = newdata3)


## ----fig.height=12, fig.width=12----------------------------------------------
oldpar <- par(no.readonly = TRUE)

par(mfrow=c(2,2))
 plot(newage,
     pred_age_mod0$ex_haz, type = "l",
     lty=1, lwd=2,
     xlab = "age",
     ylab = "excess hazard")
 
  lines(newage,
     pred_age_mod1$ex_haz, type = "l",
     lty=1, lwd=2, col = "blue")
  
legend("topleft", 
       legend = c("M0",
                  "M1"),
       lty = c(1,1),
       lwd = c(2,2),
       col = c("black", "blue"))
grid()

 plot(newage,
     pred_age_mod0$netsurv , type = "l", lty=1,
     lwd=2, xlab = "age", ylab = "net survival")
 
 lines(newage,
     pred_age_mod1$netsurv , type = "l", lty=1,
     lwd=2, col = "blue")
     grid()
legend("bottomleft", 
       legend = c("M0",
                  "M1"),
       lty = c(1,1),
       lwd = c(2,2),
       col = c("black", "blue"))
     
 plot(newage,
     pred_age_mod0$pt_cure, type = "l", lty=1, lwd=2,
     xlab = "age",
     ylab = "P(t)")
     grid()
     
 lines(newage,
     pred_age_mod1$pt_cure, type = "l", lty=1, lwd=2,
     xlab = "age", col = "blue")
     grid()
legend("bottomleft", 
       legend = c("M0",
                  "M1"),
       lty = c(1,1),
       lwd = c(2,2),
       col = c("black", "blue"))
par(oldpar)

## -----------------------------------------------------------------------------

age50 <- c(50)
agecr50 <- (age50 - mean(testiscancer$age))/sd(testiscancer$age) 

age70 <- c(70)
agecr70 <- (age70 - mean(testiscancer$age))/sd(testiscancer$age) 

time <- seq(0.1, 15, 0.1)

newdata_age50 <- with(testiscancer,
                 expand.grid(
                   event = 0,
                   age_cr = agecr50,
                   time_obs  = time))

newdata_age70 <- with(testiscancer,
                 expand.grid(
                   event = 0,
                   age_cr = agecr70,
                   time_obs  = time))


pred_age50_mod0 <- predict(object = fit_mod0,
                                  newdata = newdata_age50)

pred_age70_mod0 <- predict(object = fit_mod0,
                                  newdata = newdata_age70)



pred_age50_mod1 <- predict(object = fit_mod1,
                                  newdata = newdata_age50)

pred_age70_mod1 <- predict(object = fit_mod1,
                                  newdata = newdata_age70)



## ----fig.height=12, fig.width=12----------------------------------------------
plot(pred_age50_mod0, fun="all",conf.int=FALSE,lwd.main = 2,lwd.sub = 2)



## -----------------------------------------------------------------------------
plot(pred_age50_mod0, fun="haz",conf.int=TRUE,conf.type="plain")

plot(pred_age50_mod0, fun="surv",conf.int=TRUE,conf.type="log-log")

plot(pred_age50_mod0, fun="pt_cure",conf.int=TRUE,conf.type="log")

plot(pred_age50_mod0, fun="cumhaz",conf.int=FALSE)

plot(pred_age50_mod0, fun="logcumhaz",conf.int=TRUE,conf.type="plain")


## ----fig.height=12, fig.width=12----------------------------------------------

par(mfrow=c(2,2))
plot(
  time,
  pred_age50_mod0$ex_haz,
  type = "l",
  lty = 1,
  lwd = 2,
  xlab = "time since diagnosis (years)",
  ylab = "excess hazard", ylim = c(0,0.20))

lines(
  time,
  pred_age50_mod1$ex_haz,
  type = "l", col = "blue",
  lty = 1,
  lwd = 2)

lines(
  time,
  pred_age70_mod0$ex_haz,
  type = "l",
  lty = 2,
  lwd = 2)

lines(
  time,
  pred_age70_mod1$ex_haz,
  type = "l", col = "blue",
  lty = 2,
  lwd = 2)
grid()

legend("topright", 
       legend = c("M0 - age = 50",
                  "M1 - age = 50",
                  "M0- age = 70",
                  "M1- age = 70"),
       lty = c(1,1,2,2),
       lwd = c(2,2,2,2),
       col = c("black", "blue", "black", "blue"))
plot(time,
     pred_age50_mod0$netsurv , type = "l", lty=1,
     lwd=2,
     xlab = "time since diagnosis (years)",
     ylab = "net survival",
     ylim = c(0,1))
 
lines(time,
     pred_age50_mod1$netsurv , type = "l", lty=1,
     lwd=2, col = "blue")
     
     
     
lines(time,
     pred_age70_mod0$netsurv , type = "l", lty=2,
     lwd=2)

 
lines(time,
     pred_age70_mod1$netsurv , type = "l", lty=2,
     lwd=2, col = "blue")
     grid()
     
legend("bottomleft", 
       legend = c("M0 - age = 50",
                  "M1 - age = 50",
                  "M0- age = 70",
                  "M1- age = 70"),
       lty = c(1,1,2,2),
       lwd = c(2,2,2,2),
       col = c("black", "blue", "black", "blue"))

plot(time,
     pred_age50_mod0$pt_cure,
     type = "l", lty=1, lwd=2,
     xlab = "time since diagnosis (years)",
     ylab = "P(t)", 
     ylim = c(0,1))
     grid()
     
 lines(time,
     pred_age50_mod1$pt_cure,
     type = "l", lty=1, lwd=2,
     xlab = "age", col = "blue")

 lines(time,
     pred_age70_mod0$pt_cure,
     type = "l", lty=2, lwd=2)
 
 lines(time,
     pred_age70_mod1$pt_cure,
     type = "l", lty=2, lwd=2,
     xlab = "age", col = "blue")
     grid()
legend("bottomright", 
       legend = c("M0 - age = 50",
                  "M1 - age = 50",
                  "M0- age = 70",
                  "M1- age = 70"),
       lty = c(1,1,2,2),
       lwd = c(2,2,2,2),
       col = c("black", "blue", "black", "blue"))

par(oldpar)

## -----------------------------------------------------------------------------
date()

sessionInfo()

