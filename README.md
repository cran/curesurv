Vignette_tneh
================
Juste Goungounga, Olayidé Boussari, Laura Botta, Valérie Jooste

``` r
library(curesurv)
#> Le chargement a nécessité le package : stringr
#> Le chargement a nécessité le package : survival
```

# How to analyse survival data using Time-To-Null exccess hazard model (TNEH model)

## Additive hazard assumption in net survival setting

When the cause of death is unknown, the most common method to estimate
the cancer-related survival is net survival. Its estimation assumes that
the observed hazard $\lambda_{obs}$ is equal to the sum of the known
background mortality hazard in the general population $\lambda_{pop}$
(obtained from national Statistic Institutes such as INSEE in France)
and the excess hazard (due to cancer) $\lambda_{exc}$. For one
individual $i$, this relation can be expressed as:

$$
\lambda_{obs}(t_i|z_i) = \lambda_{pop}(t_i+a_i|z_{Pi}) + \lambda_{exc}(t_i|z_i)
$$ where $Z_P\subset Z$.

The cumulative observed hazard can be written as: $$
\Lambda_{obs}(t) = \Lambda_{pop}(t+a) + \Lambda_{exc}(t)
$$ and the net survival is obtained as following: $$
S_{n}(t) = \exp(-\Lambda_{exc}(t))
$$

# TNEH model

The TNEH model is a relatively recent excess hazard model developed by
Boussari et al. $\\$

The particularity of this model is that it enables the estimation, at
the same time as the classical parameters of a model of the excess rate,
of a quantity which is obtained by post-estimation by the classical
models: it concerns the time after which the excess rate becomes null
i.e. the cure point.

## Instantaneous excess hazard

The excess hazard proposed can be expressed as following:

$$
\lambda_{exc}(t|z;\theta) = \left(\dfrac{t}{\tau(z;\tau*)}\right)^{\alpha(z;\alpha*)-1} \left(1 - \dfrac{t}{\tau(z;\tau*)}\right)^{\beta-1} 1_{\left\{0 \le t \le \tau(z;\tau*)\right\}}
$$

where : $\\$

$\tau(z;\tau*) > 0$ is the time to cure, depends on covariates z and
vector of parameters $\tau*$. It corresponds to the vector of parameters
fitting the time-to-null excess hazard. $\\$

$\alpha(z;\alpha*) > 0$ and $\beta > 1$ are shape parameters. With
$\beta>1$, the excess hazard is forced to be null and continuous in
$\tau(z;\tau*)$. $\\$

The vector of parameters to be estimated is
$\theta = (\alpha*, \beta, \tau*)$ with $\alpha(z;\alpha*) > 0$ .

## Cumulative excess hazard

$$
\Lambda_{exc}(t|z;\theta) = \tau(z;\tau*) B \left( \alpha(z;\alpha*), \beta \right) F_{Be} \left( \dfrac{t}{\tau(z;\tau*)} ; \alpha(z;\alpha*) , \beta \right)
$$

where

B is the beta function $\\$ $F_{Be}$ is the cumulative distribution
function (cdf) of the beta distribution

## Net survival

$$
S_n(t|z) = \exp(-\Lambda_{exc}(t|z)) = \exp\left(-\tau(z;\tau*) B \left( \alpha(z;\alpha*), \beta \right) F_{Be} \left( \dfrac{t}{\tau(z;\tau*)} ; \alpha(z;\alpha*),\beta \right)\right) 
$$

## The cure fraction $\pi$

The cure fraction corresponds to the net survival at $t = \tau$ in TNEH
model. It can be expressed as:

$$
\pi(z|\theta) = \exp\left(-\Lambda_{exc}(\tau(z;\tau*)|z)\right)
        = \exp\left(-\tau(z;\tau*) B \left( \alpha(z;\alpha*), \beta \right)\right)
$$

## The probability Pi(t)

This quantity corresponds to the probability Pi(t) of being cured at a
given time t after diagnosis knowing that he/she was alive up to time t.
It can be expressed as following:

$$
Pi(t|z) = \dfrac{\pi(z|\theta)}{S_n(t|z)} 
        = \exp \left( \tau(z;\tau*) \left( B \left( \dfrac{t}{\tau(z;\tau*)} ; \alpha(z;\alpha*) , \beta \right) - B(\alpha, \beta)  \right) \right)
$$ To calculates the confidence intervals of $Pi(t|z)$, can be obtained
using the delta method. The application of this method requires the
partial derivatives of $Pi(t|z)$ with respect of the parameters of the
model. This can be written as:

$$
\dfrac{\partial Pi(t|z)}{\partial \theta} = \dfrac{1}{S_n(t|z)^2} \left( \dfrac{\partial \pi(z|\theta)}{\partial \theta} S_n(t|z) - \dfrac{\partial S_n(t|z)}{\partial \theta} \pi(z|\theta)   \right)
$$

# Fit of tneh model using R

## Without covariates

``` r
fit_ad_tneh_nocov <- curesurv(Surv(time_obs, event) ~ 1,
                             pophaz = "ehazard",
                             cumpophaz = "cumehazard",
                             model = "nmixture", dist = "tneh",  
                             link_tau = "linear",
                             data = testiscancer,
                             method_opt = "L-BFGS-B")
#> init 5 5.5 5 lower 0 1 0 upper 100 100 100 
#> next evaluation with initial values =   2
```

``` r
fit_ad_tneh_nocov
#> Call:
#> curesurv(formula = Surv(time_obs, event) ~ 1, data = testiscancer, 
#>     pophaz = "ehazard", cumpophaz = "cumehazard", model = "nmixture", 
#>     dist = "tneh", link_tau = "linear", method_opt = "L-BFGS-B")
#> 
#>          coef se(coef)      z      p
#> alpha0 2.1841   0.1032 21.166 <2e-16
#> beta   4.4413   0.5178  8.577 <2e-16
#> tau0   5.1018   0.5397  9.452 <2e-16
#> 
#> Estimates and their 95% CI after back-transformation
#>        estimates   LCI   UCI
#> alpha0     2.184 1.982 2.386
#> beta       4.441 3.426 5.456
#> tau0       5.102 4.044 6.160
#> 
#> Cured proportion exp[-ζ0* B((α0+α*Z)β)] and its 95% CI
#> 
#>    estimates    LCI    UCI
#> π0    0.8474 0.7616 0.9042
#> 
#> log-likelihood: -2633.903 (for 3 degree(s) of freedom)
#>  AIC: 5273.806
#> 
#>   n= 2000 , number of events= 949
```

``` r
newdata1 <- with(testiscancer,
  expand.grid(event = 0, time_obs  = seq(0.001,10,0.001)))
p_28 <- predict(object = fit_ad_tneh_nocov, newdata = newdata1)
```

### Plot of different estimators (hazard, survival, probability of being cured)

``` r
plot(p_28)
```

![](README_files/figure-gfm/plotpredictnocov-1.png)<!-- -->

### Cure fraction estimation precision

The confidence intervals at $1-\alpha$ level for the cure fraction $\pi$
can be written as:

$$
\left[\hat{\pi} \pm z_{1 - \alpha / 2}  \sqrt{Var(\hat{\hat{\pi}})}\right]
$$ where $$
Var(\hat{\pi}) = \dfrac{\partial \hat{\pi}}{\partial \theta} Var(\theta) \left(\dfrac{\partial \hat{\pi}}{\partial \theta}\right)^T
$$

### Time-to-cure estimation

We search the time $\text{t}=\text{TTC}_i$ from which
$\text{P}_i(t) = 1-\epsilon$. $\epsilon$ can be fixed to 0.95.

The variance formula can be expressed as:

$$
Var(TTC) = Var(g(\theta;z_i)) \simeq \left(\dfrac{\partial P(t|z_i;\theta)}{\partial t}_{|t = TTC}\right)^{-2} Var(P(t|z_i;\theta))_{|t=TTC}
$$

## With covariates acting both on parameters tau and alpha

``` r
testiscancer$age_crmin <- (testiscancer$age- min(testiscancer$age)) /sd(testiscancer$age)
```

``` r

fit_m1_ad_tneh <- curesurv(Surv(time_obs, event) ~ z_tau(age_crmin) + 
                          z_alpha(age_crmin),
                          pophaz = "ehazard",
                          cumpophaz = "cumehazard",
                          model = "nmixture", dist = "tneh",
                          link_tau = "linear",
                          data = testiscancer,
                          method_opt = "L-BFGS-B")
#> init 5 2.5 5.5 5 -4 lower 0 -5 1 0 -5 upper 100 100 100 100 100
#> Warning in diag(varcov_star): NAs introduits lors de la conversion automatique
#> Warning in diag(varcov): NAs introduits lors de la conversion automatique
#> non convergence with inititial values 1 
#> next evaluation with initial values =   2 
#> init 7.5 -1.25 7.75 2.5 3 lower 0 -5 1 0 -5 upper 100 100 100 100 100 
#> next evaluation with initial values =   3

fit_m1_ad_tneh
#> Call:
#> curesurv(formula = Surv(time_obs, event) ~ z_tau(age_crmin) + 
#>     z_alpha(age_crmin), data = testiscancer, pophaz = "ehazard", 
#>     cumpophaz = "cumehazard", model = "nmixture", dist = "tneh", 
#>     link_tau = "linear", method_opt = "L-BFGS-B")
#> 
#>                         coef se(coef)      z        p
#> alpha0               2.87737  0.24078 11.950  < 2e-16
#> alpha_1_(age_crmin) -0.50128  0.07498 -6.685 2.30e-11
#> beta                 5.15005  1.04335  4.936 7.97e-07
#> tau0                 3.25831  0.55340  5.888 3.91e-09
#> tau_1_(age_crmin)    3.46300  1.23871  2.796  0.00518
#> 
#> Estimates and their 95% CI after back-transformation
#>                     estimates   LCI   UCI
#> alpha0                  2.877 2.405 3.349
#> alpha_1_(age_crmin)     2.376 2.229 2.523
#> beta                    5.150 3.105 7.195
#> tau0                    3.258 2.174 4.343
#> tau_1_(age_crmin)       6.721 4.293 9.149
#> 
#> Cured proportion exp[-(ζ0+ζ*Z)* B((α0+α*Z)β)] and its 95% CI
#> (For each Z of (age_crmin) the others are at reference level)
#> 
#>                                estimates    LCI    UCI
#> π0z_alpha0                        0.9675 0.7989 0.9949
#> π0z_alpha(age_crmin)              0.9601 0.8183 0.9689
#> π(age_crmin)z_alpha0              0.9341 0.6231 0.9900
#> π(age_crmin)z_alpha(age_crmin)    0.9227 0.6613 0.9356
#> 
#> log-likelihood: -2544.1 (for 5 degree(s) of freedom)
#>  AIC: 5098.2
#> 
#>   n= 2000 , number of events= 949
```

``` r
  #mean of age
newdata1 <- with(testiscancer, 
                 expand.grid(event = 0,
                             age_crmin = mean(age_crmin),
                             time_obs  = seq(0.001,10,0.1)))

 pred_agemean <- predict(object = fit_m1_ad_tneh, newdata = newdata1)
```

``` r
#max of age
newdata2 <- with(testiscancer, 
                  expand.grid(event = 0,
                              age_crmin = max(age_crmin),
                              time_obs  = seq(0.001,10,0.1)))

pred_agemax <- predict(object = fit_m1_ad_tneh, newdata = newdata2)
```

``` r
# predictions at time 2 years depending on age

   newdata3 <- with(testiscancer,
      expand.grid(event = 0, 
                  age_crmin = seq(min(testiscancer$age_crmin), 
                                  max(testiscancer$age_crmin), 0.1),
                  time_obs  = 2))

pred_age_val <- predict(object = fit_m1_ad_tneh, newdata = newdata3)
```

plot of net survival for mean and maximum age

``` r


par(mfrow = c(2, 2),
    cex = 1.0)
plot(pred_agemax$time,
    pred_agemax$ex_haz,
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "Time since diagnosis",
    ylab = "excess hazard")
lines(pred_agemean$time,
     pred_agemean$ex_haz,
     type = "l",
     lty = 2,
     lwd = 2)

legend("topright",
      horiz = FALSE,
      legend = c("hE(t) age.max = 79.9", "hE(t) age.mean = 50.8"),
      col = c("black", "black"),
      lty = c(1, 2, 1, 1, 2, 2))
grid()

plot(pred_agemax$time,
    pred_agemax$netsurv,
    type = "l",
    lty = 1,
    lwd = 2,
    ylim = c(0, 1),
    xlab = "Time since diagnosis",
    ylab = "net survival")
lines(pred_agemean$time,
     pred_agemean$netsurv,
     type = "l",
     lty = 2,
     lwd = 2)
legend("bottomleft",
       horiz = FALSE,
       legend = c("Sn(t) age.max = 79.9", "Sn(t) age.mean = 50.8"),
       col = c("black", "black"),
      lty = c(1, 2, 1, 1, 2, 2))
grid()

plot(pred_agemax$time,
    pred_agemax$pt_cure,
    type = "l",
    lty = 1,
    lwd = 2,
    ylim = c(0, 1), xlim = c(0,30),
    xlab = "Time since diagnosis",
    ylab = "probability of being cured P(t)")

lines(pred_agemean$time,
     pred_agemean$pt_cure,
     type = "l",
     lty = 2,
     lwd = 2)


abline(v = pred_agemean$tau[1],
      lty = 2,
      lwd = 2,
      col = "blue")
abline(v = pred_agemean$TTC[1],
       lty = 2,
       lwd = 2,
       col = "red")
abline(v = pred_agemax$tau[1],
       lty = 1,
       lwd = 2,
       col = "blue")
abline(v = pred_agemax$TTC[1],
       lty = 1,
       lwd = 2,
      col = "red")
grid()

legend("bottomright",
       horiz = FALSE,
       legend = c("P(t) age.max = 79.9",
                 "P(t) age.mean = 50.8",
                 "TNEH age.max = 79.9",
                 "TTC age.max = 79.9",
                 "TNEH age.mean = 50.8",
                 "TTC age.mean = 50.8"),
      col = c("black", "black", "blue", "red", "blue", "red"),
      lty = c(1, 2, 1, 1, 2, 2))
```

![](README_files/figure-gfm/predict_tnehcov-1.png)<!-- -->

``` r
val_age <- seq(min(testiscancer$age_crmin),
               max(testiscancer$age_crmin),
               0.1) * sd(testiscancer$age) +  min(testiscancer$age)


pred_age_val <- predict(object = fit_m1_ad_tneh, newdata = newdata3)
```

``` r
par(mfrow=c(2,2))
 plot(val_age,
     pred_age_val$ex_haz, type = "l",
     lty=1, lwd=2,
     xlab = "age",
     ylab = "excess hazard")
grid()

 plot(val_age,
     pred_age_val$netsurv, type = "l", lty=1,
     lwd=2, xlab = "age", ylab = "net survival")
     grid()

 plot(val_age,
     pred_age_val$pt_cure, type = "l", lty=1, lwd=2,
     xlab = "age",
     ylab = "Pi(t)")
     grid()
par(mfrow=c(1,1))
```

![](README_files/figure-gfm/predicted_to_plot_by_age-1.png)<!-- -->

## With covariates acting only parameters adjusting the parameter of the time-to-null excess hazard tau only

``` r

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
#> init 5 5.5 5 -4 lower 0 1 0 -5 upper 100 100 100 100
#> Warning in diag(varcov_star): NAs introduits lors de la conversion automatique
#> Warning in diag(varcov): NAs introduits lors de la conversion automatique
#> non convergence with inititial values 1 
#> next evaluation with initial values =   2 
#> init 7.5 3.25 7.5 -11 lower 0 1 0 -5 upper 100 100 100 100
#> Warning in diag(varcov_star): NAs introduits lors de la conversion automatique

#> Warning in diag(varcov_star): NAs introduits lors de la conversion automatique
#> non convergence with inititial values 2 
#> next evaluation with initial values =   3 
#> init 2.5 7.75 2.5 3 lower 0 1 0 -5 upper 100 100 100 100
#> Warning in diag(varcov_star): NAs introduits lors de la conversion automatique

#> Warning in diag(varcov_star): NAs introduits lors de la conversion automatique
#> non convergence with inititial values 3 
#> next evaluation with initial values =   4 
#> init 3.75 4.375 6.25 -14.5 lower 0 1 0 -5 upper 100 100 100 100
#> Warning in diag(varcov_star): NAs introduits lors de la conversion automatique

#> Warning in diag(varcov_star): NAs introduits lors de la conversion automatique
#> non convergence with inititial values 4 
#> next evaluation with initial values =   5 
#> init 8.75 8.875 1.25 -0.5 lower 0 1 0 -5 upper 100 100 100 100 
#> next evaluation with initial values =   6
```

``` r
fit_ad_tneh_covtau
#> Call:
#> curesurv(formula = Surv(time_obs, event) ~ z_tau(age_cr), data = testiscancer, 
#>     pophaz = "ehazard", cumpophaz = "cumehazard", model = "nmixture", 
#>     dist = "tneh", link_tau = "linear", method_opt = "L-BFGS-B")
#> 
#>                  coef se(coef)      z        p
#> alpha0         1.9753   0.1299 15.206  < 2e-16
#> beta           5.3066   1.1286  4.702 2.58e-06
#> tau0           7.4380   1.8109  4.107 4.00e-05
#> tau_1_(age_cr) 2.3159   0.8529  2.715  0.00662
#> 
#> Estimates and their 95% CI after back-transformation
#>                estimates   LCI    UCI
#> alpha0             1.975 1.721  2.230
#> beta               5.307 3.095  7.519
#> tau0               7.438 3.889 10.987
#> tau_1_(age_cr)     9.754 8.082 11.426
#> 
#> Cured proportion exp[-(ζ0+ζ*Z)* B((α0+α*Z)β)] and its 95% CI
#> (For each Z of (age_cr) the others are at reference level)
#> 
#>           Estimates    LCI    UCI
#> π0            0.794 0.6535 0.8909
#> π(age_cr)     0.739 0.4130 0.8868
#> 
#> log-likelihood: -2610.768 (for 4 degree(s) of freedom)
#>  AIC: 5229.537
#> 
#>   n= 2000 , number of events= 949
```

``` r
newdata2 <- with(testiscancer,
                 expand.grid(event = 0, 
                             time_obs  = seq(0.001, 10, 0.001),
                             age_cr = c(-0.9577, -0.2751, 0.2849) ))
newdata2_1stqu <- newdata2[newdata2$age_cr==-0.9577,]
newdata2_2rdqu <- newdata2[newdata2$age_cr==-0.2751,]
newdata2_3rdqu <- newdata2[newdata2$age_cr==0.2849,]

p1stqu <- predict(object = fit_ad_tneh_covtau, newdata = newdata2_1stqu)
p2rdqu <- predict(object = fit_ad_tneh_covtau, newdata = newdata2_2rdqu)
p3rdqu <- predict(object = fit_ad_tneh_covtau, newdata = newdata2_3rdqu)
```

``` r
oldpar <- par(no.readonly = FALSE)
par(mfrow = c(2,2))
plot(p1stqu, 
     main = "Excess hazard for age 20", 
     fun = "haz")
plot(p2rdqu,
     fun = "haz",
     main = "Excess hazard for age 51")
plot(p3rdqu,
     fun = "haz",
     main = "Excess hazard for age 69")
par(mfrow = c(1,1))
```

![](README_files/figure-gfm/plot_for_specific_age-1.png)<!-- -->

``` r
par(oldpar)
#> Warning in par(oldpar): le paramètre graphique "cin" ne peut être changé
#> Warning in par(oldpar): le paramètre graphique "cra" ne peut être changé
#> Warning in par(oldpar): le paramètre graphique "csi" ne peut être changé
#> Warning in par(oldpar): le paramètre graphique "cxy" ne peut être changé
#> Warning in par(oldpar): le paramètre graphique "din" ne peut être changé
#> Warning in par(oldpar): le paramètre graphique "page" ne peut être changé
```

## With covariates acting only on scale parameter alpha

``` r

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
#> init 5 2.5 5.5 5 lower 0 -5 1 0 upper 100 100 100 100 
#> next evaluation with initial values =   2
```

``` r
fit_ad_tneh_covalpha
#> Call:
#> curesurv(formula = Surv(time_obs, event) ~ z_alpha(age_cr), data = testiscancer, 
#>     pophaz = "ehazard", cumpophaz = "cumehazard", model = "nmixture", 
#>     dist = "tneh", link_tau = "linear", method_opt = "L-BFGS-B")
#> 
#>                      coef se(coef)      z        p
#> alpha0            2.06862  0.11152 18.550  < 2e-16
#> alpha_1_(age_cr) -0.46785  0.06331 -7.389 1.48e-13
#> beta              4.77703  0.81573  5.856 4.74e-09
#> tau0              6.09881  1.11797  5.455 4.89e-08
#> 
#> Estimates and their 95% CI after back-transformation
#>                  estimates   LCI   UCI
#> alpha0               2.069 1.850 2.287
#> alpha_1_(age_cr)     1.601 1.477 1.725
#> beta                 4.777 3.178 6.376
#> tau0                 6.099 3.908 8.290
#> 
#> Cured proportion exp[-ζ0* B((α0+α*Z)β)] and its 95% CI
#> 
#>           estimates    LCI    UCI
#> π0           0.8181 0.4760 0.9485
#> π(age_cr)    0.7709 0.4124 0.7536
#> 
#> log-likelihood: -2586.138 (for 4 degree(s) of freedom)
#>  AIC: 5180.275
#> 
#>   n= 2000 , number of events= 949
```

``` r
p4_28 <- predict(object = fit_ad_tneh_covalpha,
                 newdata = newdata2_1stqu)
p4_50 <- predict(object = fit_ad_tneh_covalpha,
                 newdata = newdata2_2rdqu)
p4_74 <- predict(object = fit_ad_tneh_covalpha,
                 newdata = newdata2_2rdqu)
```

### plot of estimation of probability Pi(t) of being cured at a given time t after diagnosis knowing that he/she was alive up to time t

``` r

#| echo: true
#| message: false
#| warning: false
#| include: true
#| fig.height: 15
#| fig.width: 15
par(mfrow = c(2,2))
plot(p4_28, fun = "pt_cure")
plot(p4_50, fun = "pt_cure")
plot(p4_74, fun = "pt_cure")
par(mfrow = c(1,1))
```

![](README_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->
