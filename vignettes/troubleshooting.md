---
title: "Troubleshooting with glmmTMB"
date: "2018-06-20"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{troubleshooting}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



This vignette covers common problems that occur while using `glmmTMB`. 
The contents will expand with experience.

If your problem is not covered below, try updating to the latest version of `glmmTMB` on GitHub. The developers might have solved the problem in a newer version. 

#Warnings

## Model convergence problem; non-positive-definite Hessian matrix

You may see the same warning as in the following example:

```r
zinbm0 = glmmTMB(count~spp + (1|site), zi=~spp, Salamanders, family=nbinom2)
```

```
## Warning in fitTMB(TMBStruc): Model convergence problem; non-positive-
## definite Hessian matrix. See vignette('troubleshooting')
```

This error states that the point that `glmmTMB` has identified as the putative maximum-likelihood estimate, the curvature of the log-likelihood surface does not seem to be consistent with `glmmTMB` really having found a maximum: instead, it is upward-curving, or flat, in some direction(s).

It can occur:

- when a model is overparameterized (i.e. the data does not contain information to estimate the parameters)
- when a random-effect variance is estimated to be zero, or random-effect terms are estimated to be perfectly correlated (often caused by having too few levels of the random-effect grouping variable)
- when *complete separation* occurs in a binomial model: some categories in the model contain proportions that are either all 0 or all 1

How do we diagnose the problem?
First, look at the estimated coefficients and see if any of them are extreme. If you're using a non-identity
link function (e.g. log, logit), then parameter values with $|\beta|>10$ are suspect (for logit link, this
implies probabilities very close to 0 or 1; for log link, this implies mean counts that are
close to 0 or gigantic).

Inspecting the fixed-effect estimates for this model:


```r
fixef(zinbm0)
```

```
## 
## Conditional model:
## (Intercept)        sppPR        sppDM      sppEC-A      sppEC-L  
##     -0.5377      -0.6531       0.3303      -0.2001       0.6615  
##    sppDES-L        sppDF  
##      0.7993       0.3714  
## 
## Zero-inflation model:
## (Intercept)        sppPR        sppDM      sppEC-A      sppEC-L  
##     -17.825       17.852       -5.202       17.486       15.431  
##    sppDES-L        sppDF  
##      14.708       15.485
```

The zero-inflation intercept parameter is tiny (-17): since the parameters
are estimated on the logit scale, we back-transform with `plogis(-17)` to see the at the zero-inflation probability for the baseline level is about $4 \times 10^{-8}$)). Many of the other ZI parameters are very large, compensating for the intercept: the estimated zero-inflation probabilities for all species are


```r
ff <- fixef(zinbm0)$zi
round(plogis(c(sppGP=unname(ff[1]),ff[-1]+ff[1])),3)
```

```
##    sppGP    sppPR    sppDM  sppEC-A  sppEC-L sppDES-L    sppDF 
##    0.000    0.507    0.000    0.416    0.084    0.042    0.088
```

Since the baseline probability is already effectively zero,
making the intercept parameter larger or smaller will have very little effect - the likelihood is flat,
which leads to the non-positive-definite warning.

Now that we suspect the problem is in the zero-inflation component,
we can try to come up with ways of simplifying the model:
for example, we could use a model that compared the first species ("GP") to the rest:


```r
Salamanders <- transform(Salamanders, GP=as.numeric(spp=="GP"))
zinbm0_A = update(zinbm0, ziformula=~GP)
```

This fits without a warning, although the GP zero-inflation parameter is still extreme:


```r
fixef(zinbm0_A)[["zi"]]
```

```
## (Intercept)          GP 
##   -2.890515  -15.414542
```

Another possibility would be to fit the variation among species in the zero-inflation parameter
as a random effect, rather than a fixed effect: this is slightly more parsimonious.
This again fits without an error, although both the average level of
zero-inflation and the among-species variation are estimated as very small:


```r
zinbm0_B = update(zinbm0, ziformula=~(1|spp))
fixef(zinbm0_B)[["zi"]]
```

```
## (Intercept) 
##   -16.54352
```

```r
VarCorr(zinbm0_B)
```

```
## 
## Conditional model:
##  Groups Name        Std.Dev.
##  site   (Intercept) 1.3894  
## 
## Zero-inflation model:
##  Groups Name        Std.Dev. 
##  spp    (Intercept) 0.0078215
```

<!-- FIXME: updating here does weird things
zinbm1 = update(zinbm0, ziformula=~mined, Salamanders, family=nbinom2)
-->

The original analysis considered variation in zero-inflation by site status
(mined or not mined) rather than by species - this simpler model only tries
to estimate two parameters (mined + difference between mined and no-mining)
rather than 7 (one per species) for the zero-inflation model.


```r
zinbm1 = glmmTMB(count~spp + (1|site), zi=~mined, Salamanders, family=nbinom2)
fixef(zinbm1)[["zi"]]
```

```
## (Intercept)     minedno 
##   0.3787986 -17.5118413
```

This again fits without a warning, but we see that the zero-inflation is effectively
zero in the unmined ("minedno") condition (`plogis(0.38-17.5)` is
approximately $$4 \times 10^{-8}$$). We can estimate the confidence interval, but
it takes some extra work: the default Wald standard errors and confidence intervals
are useless in this case.


```r
## at present we need to specify the parameter by number; for
##  extreme cases need to specify the parameter range
## (not sure why the upper bound needs to be so high ... ?)
cc = confint(zinbm1,method="uniroot",parm=9, parm.range=c(-20,20))
print(cc)
```

```
##            2.5 %    97.5 %  Estimate
## zi~minedno    NA -2.083725 -17.51184
```

The lower CI is not defined; the upper CI is -2.08, i.e. we can state
that the zero-inflation probability is less than `plogis(-2.08)` = 0.11.

More broadly, general inspection of the data (e.g., plotting the response against potential covariates)
should help to diagnose overly complex models.

In some cases, scaling predictor variables may help.

In general models with non-positive definite Hessian matrices should be excluded from further consideration.

##Model convergence problem:  eigenvalue problems


```r
m1 = glmmTMB(count~spp + mined + (1|site), zi=~spp + mined, Salamanders, family=genpois)
```

```
## Warning in fitTMB(TMBStruc): Model convergence problem; extreme or very
## small eigen values detected. See vignette('troubleshooting')
```

In this example, the fixed-effect covariance matrix is `NaN`. It may have to do with the generalized Poisson (`genpois`) distribution, which is known to have convergence problems; luckily, the negative binomial (`nbinom1` and `nbinom2`) and/or Conway-Maxwell Poisson (`compois`) are good alternatives. 

Models with convergence problems should be excluded from further consideration, in general.

In some cases, extreme eigenvalues may be caused by having predictor variables that are on very different scales: try rescaling, and centering, continuous predictors in the model.

##NA/NaN function evaluation

> Warning in nlminb(start = par, objective = fn, gradient = gr) : NA/NaN function evaluation

This warning occurs when the optimizer visits a region of parameter space that is invalid. It is not a problem as long as the optimizer has left that region of parameter space upon convergence, which is indicated by an absence of the model convergence warnings described above. 

The following warnings indicate possibly-transient numerical problems with the fit, and can be treated in the same way (i.e. ignored if there are no errors or convergence warnings about the final fitted model).


```r
Cholmod warning 'matrix not positive definite'
```

```r
Warning in f(par, order = order, ...) : value out of range in 'lgamma'
```
# Errors

## NA/NaN gradient evaluation


```r
dat1 = expand.grid(y=-1:1, rep=1:10)
m1 = glmmTMB(y~1, dat1, family=nbinom2)
```

```
## Error in nlminb(start = par, objective = fn, gradient = gr, control = control$optCtrl): NA/NaN gradient evaluation
```

```
## Timing stopped at: 0.002 0 0.002
```
The error occurs here because the negative binomial distribution is not appropriate for data with negative values.

If you see this error, check that the response variable meets the assumptions of the specified distribution.

## gradient length


> Error in nlminb(start = par, objective = fn, gradient = gr) : gradient function must return a numeric vector of length x

> Error in optimHess(par.fixed, obj$fn, obj$gr): gradient in optim evaluated to length x

Try rescaling predictor variables. Try a simpler model and build up.


