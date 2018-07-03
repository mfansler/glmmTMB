---
title: "Covariance structures with glmmTMB"
author: "Kasper Kristensen"
date: "2018-06-19"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{random effect structures}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



This vignette demonstrates some of the covariance structures available in the `glmmTMB` package.
Currently the available covariance structures are:

| Covariance                       | Notation      | Parameter count | Requirement |
|----------------------------------|---------------|-----------------|-------------|
| Heterogeneous unstructured       | `us`          |  $n(n+1)/2$     |             |
| Heterogeneous Toeplitz           | `toep`        |  $2n-1$         |             |
| Heterogeneous compound symmetry  | `cs`          |  $n+1$          |             |
| Heterogeneous diagonal           | `diag`        |  $n$            |             |
| AR(1)                            | `ar1`         |  $2$            |             |
| Ornstein–Uhlenbeck               | `ou`          |  $2$            | Coordinates |
| Spatial exponential              | `exp`         |  $2$            | Coordinates |
| Spatial Gaussian                 | `gau`         |  $2$            | Coordinates |
| Spatial Matern                   | `mat`         |  $3$            | Coordinates |

The word 'heterogeneous' refers to the marginal variances of the
model. Beyond correlation parameters, a heterogeneous structure uses
$n$ additional variance parameters where $n$ is the dimension.

Some of the structures require temporal or spatial coordinates. We
will show examples of this in a later section.

## The AR(1) covariance structure

### Demonstration on simulated data

First, let's consider a simple time series model. Assume that our
measurements $Y(t)$ are given at discrete times $t \in \{1,...,n\}$ by

$$Y(t) = \mu + X(t) + \varepsilon(t)$$

where

- $\mu$ is the mean value parameter.
- $X(t)$ is a stationary AR(1) process, i.e. has covariance $cov(X(s),
  X(t)) = \sigma^2\exp(-\theta |t-s|)$.
- $\varepsilon(t)$ is iid. $N(0,\sigma_0^2)$ measurement error.

A simulation experiment is set up using the parameters

| Description            | Parameter     | Value |
|------------------------|---------------|-------|
| Mean                   | $\mu$         | 0     |
| Process variance       | $\sigma^2$    | 1     |
| Measurement variance   | $\sigma_0^2$  | 1     |
| One-step correlation   | $e^{-\theta}$ | 0.7   |

The following R-code draws a simulation based on these parameter
values.  For illustration purposes we consider a very short time
series.


```r
n <- 6                                              ## Number of time points
x <- mvrnorm(mu = rep(0,n),
             Sigma = .7 ^ as.matrix(dist(1:n)) )    ## Simulate the process using the MASS package
y <- x + rnorm(n)                                   ## Add measurement noise
```

In order to fit the model with `glmmTMB` we must first specify a time
variable as a *factor*. The factor *levels* correspond to unit spaced
time points.


```r
times <- factor(1:n)
levels(times)
```

```
## [1] "1" "2" "3" "4" "5" "6"
```

We also need a grouping variable. In the current case there is only
one time-series so the grouping is:


```r
group <- factor(rep(1,n))
```

We combine the data into a single data frame (not absolutely
required, but good practice):


```r
dat0 <- data.frame(y,times,group)
```

Now fit the model using


```r
glmmTMB(y ~ ar1(times + 0 | group), data=dat0)
```

This formula notation follows that of the `lme4` package.

- The left hand side of the bar `times + 0` corresponds to a design
  matrix $Z$ linking observation vector $y$ (rows) with a random
  effects vector $u$ (columns).
- The distribution of $u$ is `ar1` (this is the only `glmmTMB`
  specific part of the formula).
- The right hand side of the bar splits the above specification
  independently among groups. Each group has its own separate $u$
  vector but shares the same parameters for the covariance structure.

After running the model, we find the parameter estimates $\mu$
(intercept), $\sigma_0^2$ (dispersion), $\sigma$ (Std. Dev.) and
$e^{-\theta}$ (First off-diagonal of "Corr") in the output:

> FIXME: Try a longer time series when the print.VarCorr is fixed.


```
## Formula:          y ~ ar1(times + 0 | group)
## Data: dat0
##      AIC      BIC   logLik df.resid 
## 18.76343 17.93047 -5.38171        2 
## Random-effects (co)variances:
## 
## Conditional model:
##  Groups   Name   Std.Dev.  Corr       
##  group    times1 7.018e-01 -0.58 (ar1)
##  Residual        1.165e-05            
## 
## Number of obs: 6 / Conditional model: group, 1
## 
## Dispersion estimate for gaussian family (sigma^2): 1.36e-10 
## 
## Fixed Effects:
## 
## Conditional model:
## (Intercept)  
##       1.072
```

### Increasing the sample size

A single time series of 6 time points is not sufficient to identify
the parameters. We could either increase the length of the time series
or increase the number of groups. We'll try the latter:


```r
simGroup <- function(g) {
    x <- mvrnorm(mu = rep(0,n),
             Sigma = .7 ^ as.matrix(dist(1:n)) )    ## Simulate the process
    y <- x + rnorm(n)                               ## Add measurement noise
    times <- factor(1:n)
    group <- factor(rep(g,n))
    data.frame(y, times, group)
}
simGroup(1)
```

```
##           y times group
## 1 0.3008231     1     1
## 2 0.1835385     2     1
## 3 0.3750261     3     1
## 4 1.5819573     4     1
## 5 1.5776466     5     1
## 6 0.2072020     6     1
```

Generate a dataset with 1000 groups:


```r
dat1 <- do.call("rbind", lapply(1:1000, simGroup) )
```

And fitting the model on this larger dataset gives estimates close to
the true values:


```r
(fit.ar1 <- glmmTMB(y ~ ar1(times + 0 | group), data=dat1))
```

```
## Formula:          y ~ ar1(times + 0 | group)
## Data: dat1
##       AIC       BIC    logLik  df.resid 
##  20549.90  20576.70 -10270.95      5996 
## Random-effects (co)variances:
## 
## Conditional model:
##  Groups   Name   Std.Dev. Corr      
##  group    times1 1.003    0.72 (ar1)
##  Residual        1.019              
## 
## Number of obs: 6000 / Conditional model: group, 1000
## 
## Dispersion estimate for gaussian family (sigma^2): 1.04 
## 
## Fixed Effects:
## 
## Conditional model:
## (Intercept)  
##     0.04378
```

## The unstructured covariance

We can try to fit an unstructured covariance to the previous dataset
`dat`. For this case an unstructured covariance has 15
correlation parameters and 6 variance parameters. Adding
$\sigma_0^2 I$ on top would cause a strict
overparameterization. Hence, when fitting the model with `glmmTMB`, we
have to disable the $\varepsilon$ term (the dispersion):


```r
fit.us <- glmmTMB(y ~ us(times + 0 | group), data=dat1, dispformula=~0)
fit.us$sdr$pdHess ## Converged ?
```

```
## [1] TRUE
```

The estimated variance and correlation parameters are:


```r
VarCorr(fit.us)
```

```
## 
## Conditional model:
##  Groups   Name   Std.Dev.   Corr                          
##  group    times1 1.41313653                               
##           times2 1.38976885 0.353                         
##           times3 1.50181088 0.273 0.384                   
##           times4 1.40399008 0.166 0.216 0.352             
##           times5 1.43307431 0.160 0.192 0.239 0.354       
##           times6 1.43753488 0.108 0.136 0.195 0.273 0.345 
##  Residual        0.00012207
```

The estimated correlation is approximately constant along diagonals
(apparent Toeplitz structure) and we note that the first off-diagonal
is now ca. half the true value (0.7) because the dispersion is
effectively included in the estimated covariance matrix.

## The Toeplitz structure

The next natural step would be to reduce the number of parameters by
collecting correlation parameters within the same off-diagonal. This
amounts to 5 correlation parameters and 6 variance
parameters.

> FIXME: Explain why dispformula=~1 causes over-parameterization


```r
fit.toep <- glmmTMB(y ~ toep(times + 0 | group), data=dat1, dispformula=~0)
fit.toep$sdr$pdHess ## Converged ?
```

```
## [1] TRUE
```

The estimated variance and correlation parameters are:


```r
VarCorr(fit.toep)
```

```
## 
## Conditional model:
##  Groups   Name   Std.Dev.   Corr                          
##  group    times1 1.41255279                               
##           times2 1.38784528 0.358                         
##           times3 1.49610367 0.249 0.358                   
##           times4 1.40698827 0.188 0.249 0.358             
##           times5 1.43694789 0.148 0.188 0.249 0.358       
##           times6 1.43875115 0.106 0.148 0.188 0.249 0.358 
##  Residual        0.00012207
```

The residual variance appears downward biased. REML estimation
(currently not part of `glmmTMB`) would probably give a better
estimate of the variance and thereby the correlation parameters.

> FIXME: Add REML argument to glmmTMB

## Compound symmetry

The compound symmetry structure collects all off-diagonal elements of
the correlation matrix to one common value.

> FIXME: Explain why dispformula=~1 causes over-parameterization


```r
fit.cs <- glmmTMB(y ~ cs(times + 0 | group), data=dat1, dispformula=~0)
fit.cs$sdr$pdHess ## Converged ?
```

```
## [1] TRUE
```

The estimated variance and correlation parameters are:


```r
VarCorr(fit.cs)
```

```
## 
## Conditional model:
##  Groups   Name   Std.Dev.   Corr      
##  group    times1 1.43122953 0.250 (cs)
##  Residual        0.00012207
```

## Anova tables

The models ar1, toep, and us are nested so we can use:


```r
anova(fit.ar1, fit.toep, fit.us)
```

```
## Data: dat1
## Models:
## fit.ar1: y ~ ar1(times + 0 | group), zi=~0, disp=~1
## fit.toep: y ~ toep(times + 0 | group), zi=~0, disp=~0
## fit.us: y ~ us(times + 0 | group), zi=~0, disp=~0
##          Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
## fit.ar1   4 20550 20577 -10271    20542                         
## fit.toep 12 20556 20637 -10266    20532 9.4895      8     0.3027
## fit.us   22 20570 20717 -10263    20526 6.5103     10     0.7707
```

The model `cs` is a sub-model of `toep`:


```r
anova(fit.cs, fit.toep)
```

```
## Data: dat1
## Models:
## fit.cs: y ~ cs(times + 0 | group), zi=~0, disp=~0
## fit.toep: y ~ toep(times + 0 | group), zi=~0, disp=~0
##          Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
## fit.cs    8 20703 20756 -10343    20687                             
## fit.toep 12 20556 20637 -10266    20532 154.48      4  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Adding coordinate information

Coordinate information can be added to a variable using the `glmmTMB`
function `numFactor`. This is necessary in order to use those
covariance structures that require coordinates. For example, if we
have the numeric coordinates


```r
x <- sample(1:2, 10, replace=TRUE)
y <- sample(1:2, 10, replace=TRUE)
```

we can generate a factor representing $(x,y)$ coordinates by


```r
(pos <- numFactor(x,y))
```

```
##  [1] (2,1) (2,2) (2,1) (1,2) (2,1) (1,1) (1,2) (2,2) (1,1) (1,2)
## Levels: (1,1) (2,1) (1,2) (2,2)
```

Numeric coordinates can be recovered from the factor levels:


```r
parseNumLevels(levels(pos))
```

```
##      [,1] [,2]
## [1,]    1    1
## [2,]    2    1
## [3,]    1    2
## [4,]    2    2
```

In order to try the remaining structures on our test data we
re-interpret the time factor using `numFactor`:


```r
dat1$times <- numFactor(dat1$times)
levels(dat1$times)
```

```
## [1] "(1)" "(2)" "(3)" "(4)" "(5)" "(6)"
```

## Ornstein–Uhlenbeck

Having the numeric times encoded in the factor levels we can now try
the Ornstein–Uhlenbeck covariance structure.


```r
fit.ou <- glmmTMB(y ~ ou(times + 0 | group), data=dat1)
fit.ou$sdr$pdHess ## Converged ?
```

```
## [1] TRUE
```

It should give the exact same results as `ar1` in this case since the
times are equidistant:


```r
VarCorr(fit.ou)
```

```
## 
## Conditional model:
##  Groups   Name     Std.Dev. Corr                          
##  group    times(1) 1.0034                                 
##           times(2) 1.0034   0.721                         
##           times(3) 1.0034   0.520 0.721                   
##           times(4) 1.0034   0.375 0.520 0.721             
##           times(5) 1.0034   0.270 0.375 0.520 0.721       
##           times(6) 1.0034   0.195 0.270 0.375 0.520 0.721 
##  Residual          1.0187
```

However, note the differences between `ou` and `ar1`:

- `ou` can handle irregular time points.
- `ou` only allows positive correlation between neighboring time points.

## Spatial correlations

The structures `exp`, `gau` and `mat` are meant to used for spatial
data. They all require a Euclidean distance matrix which is calculated
internally based on the coordinates. Here, we will try these models on
the simulated time series data.

An example with spatial data is presented in a later section.

### Matern


```r
fit.mat <- glmmTMB(y ~ mat(times + 0 | group), data=dat1, dispformula=~0)
fit.mat$sdr$pdHess ## Converged ?
```

```
## [1] TRUE
```


```r
VarCorr(fit.mat)
```

```
## 
## Conditional model:
##  Groups   Name     Std.Dev.   Corr                          
##  group    times(1) 1.42996896                               
##           times(2) 1.42996896 0.357                         
##           times(3) 1.42996896 0.250 0.357                   
##           times(4) 1.42996896 0.187 0.250 0.357             
##           times(5) 1.42996896 0.143 0.187 0.250 0.357       
##           times(6) 1.42996896 0.111 0.143 0.187 0.250 0.357 
##  Residual          0.00012207
```

### Gaussian

"Gaussian" refers here to a Gaussian decay in correlation with distance,
i.e. $\rho = \exp(-d x^2)$, not to the conditional distribution ("family").


```r
fit.gau <- glmmTMB(y ~ gau(times + 0 | group), data=dat1, dispformula=~0)
fit.gau$sdr$pdHess ## Converged ?
```

```
## [1] TRUE
```


```r
VarCorr(fit.gau)
```

```
## 
## Conditional model:
##  Groups   Name     Std.Dev.   Corr                          
##  group    times(1) 1.41749843                               
##           times(2) 1.41749843 0.280                         
##           times(3) 1.41749843 0.006 0.280                   
##           times(4) 1.41749843 0.000 0.006 0.280             
##           times(5) 1.41749843 0.000 0.000 0.006 0.280       
##           times(6) 1.41749843 0.000 0.000 0.000 0.006 0.280 
##  Residual          0.00012207
```

### Exponential


```r
fit.exp <- glmmTMB(y ~ exp(times + 0 | group), data=dat1)
fit.exp$sdr$pdHess ## Converged ?
```

```
## [1] TRUE
```


```r
VarCorr(fit.exp)
```

```
## 
## Conditional model:
##  Groups   Name     Std.Dev. Corr                          
##  group    times(1) 1.0034                                 
##           times(2) 1.0034   0.721                         
##           times(3) 1.0034   0.520 0.721                   
##           times(4) 1.0034   0.375 0.520 0.721             
##           times(5) 1.0034   0.270 0.375 0.520 0.721       
##           times(6) 1.0034   0.195 0.270 0.375 0.520 0.721 
##  Residual          1.0187
```

### A spatial covariance example

Starting out with the built in `volcano` dataset we reshape it to a
`data.frame` with pixel intensity `z` and pixel position `x` and `y`:


```r
d <- data.frame(z = as.vector(volcano),
                x = as.vector(row(volcano)),
                y = as.vector(col(volcano)))
```

Next, add random normal noise to the pixel intensities and extract a
small subset of 100 pixels. This is our spatial dataset:


```r
set.seed(1)
d$z <- d$z + rnorm(length(volcano), sd=15)
d <- d[sample(nrow(d), 100), ]
```

Display sampled noisy volcano data:


```r
volcano.data <- array(NA, dim(volcano))
volcano.data[cbind(d$x, d$y)] <- d$z
image(volcano.data, main="Spatial data")
```

![plot of chunk volcano_data_image](figure/volcano_data_image-1.png)

Based on this data, we'll attempt to re-construct the original image.

As model, it is assumed that the original image `image(volcano)` is a
realization of a random field with correlation decaying exponentially
with distance between pixels.

Denoting by $u(x,y)$ this random field the model for the observations is

\[ z_{i} = \mu + u(x_i,y_i) + \varepsilon_i \]

To fit the model, a `numFactor` and a dummy grouping variable must be
added to the dataset:


```r
d$pos <- numFactor(d$x, d$y)
d$group <- factor(rep(1, nrow(d)))
```

The model is fit by


```r
f <- glmmTMB(z ~ 1 + exp(pos + 0 | group), data=d)
```

Recall that a standard deviation `sd=15` was used to distort the
image. A confidence interval for this parameter is


```r
confint(f, "sigma")
```

```
##          2.5 %   97.5 % Estimate
## sigma 9.492257 16.46098 12.50007
```

The glmmTMB `predict` method can predict unseen levels of the random
effects. For instance to predict a 3-by-3 corner of the image one
could construct the new data:


```r
newdata <- data.frame( pos=numFactor(expand.grid(x=1:3,y=1:3)) )
newdata$group <- factor(rep(1, nrow(newdata)))
newdata
```

```
##     pos group
## 1 (1,1)     1
## 2 (2,1)     1
## 3 (3,1)     1
## 4 (1,2)     1
## 5 (2,2)     1
## 6 (3,2)     1
## 7 (1,3)     1
## 8 (2,3)     1
## 9 (3,3)     1
```

and predict using


```r
predict(f, newdata, type="response", allow.new.levels=TRUE)
```

```
## [1] 104.4634 104.8042 105.1246 104.6228 104.9824 105.3195 104.7828 105.1666
## [9] 105.5366
```

A specific image column can thus be predicted using the function


```r
predict_col <- function(i) {
    newdata <- data.frame( pos = numFactor(expand.grid(1:87,i)))
    newdata$group <- factor(rep(1,nrow(newdata)))
    predict(f, newdata=newdata, type="response", allow.new.levels=TRUE)
}
```

Prediction of the entire image is carried out by (this takes a while...):


```r
pred <- sapply(1:61, predict_col)
```

Finally plot the re-constructed image by


```r
image(pred, main="Reconstruction")
```

![plot of chunk image_results](figure/image_results-1.png)

## Mappings

For various advanced purposes, such as computing likelihood profiles, it is useful
to know the details of the parameterization of the models - the scale on which
the parameters are defined (e.g. standard deviation, variance, or log-standard deviation
for variance parameters) and their order.

### Unstructured

For an unstructured matrix of size `n` parameters `1:n` represent the log-standard deviations while the remaining `n(n-1)/2` (i.e. `(n+1):(n:(n*(n+1)/2))`) are the elements of the Cholesky factor of the correlation matrix, filled in row-wise order (see [TMB documentation](http://kaskr.github.io/adcomp/classUNSTRUCTURED__CORR__t.html))


```r
vv0 <- VarCorr(fit.us)
vv1 <- vv0$cond$group          ## extract 'naked' V-C matrix
n <- nrow(vv1)
rpars <- getME(fit.us,"theta") ## extract V-C parameters
## first n parameters are log-std devs:
all.equal(unname(diag(vv1)),exp(rpars[1:n])^2)
```

```
## [1] TRUE
```

```r
## now try correlation parameters:
cpars <- rpars[-(1:n)]
length(cpars)==n*(n-1)/2      ## the expected number
```

```
## [1] TRUE
```

```r
cc <- diag(n)
cc[upper.tri(cc)] <- cpars
L <- crossprod(cc)
D <- diag(1/sqrt(diag(L)))
D %*% L %*% D
```

```
##           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
## [1,] 1.0000000 0.3534513 0.2727986 0.1655881 0.1600080 0.1076348
## [2,] 0.3534513 1.0000000 0.3840508 0.2164539 0.1920130 0.1363969
## [3,] 0.2727986 0.3840508 1.0000000 0.3519413 0.2394706 0.1949200
## [4,] 0.1655881 0.2164539 0.3519413 1.0000000 0.3543410 0.2731236
## [5,] 0.1600080 0.1920130 0.2394706 0.3543410 1.0000000 0.3445279
## [6,] 0.1076348 0.1363969 0.1949200 0.2731236 0.3445279 1.0000000
```

```r
unname(attr(vv1,"correlation"))
```

```
##           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
## [1,] 1.0000000 0.3534513 0.2727986 0.1655881 0.1600080 0.1076348
## [2,] 0.3534513 1.0000000 0.3840508 0.2164539 0.1920130 0.1363969
## [3,] 0.2727986 0.3840508 1.0000000 0.3519413 0.2394706 0.1949200
## [4,] 0.1655881 0.2164539 0.3519413 1.0000000 0.3543410 0.2731236
## [5,] 0.1600080 0.1920130 0.2394706 0.3543410 1.0000000 0.3445279
## [6,] 0.1076348 0.1363969 0.1949200 0.2731236 0.3445279 1.0000000
```

> FIXME: why are these not quite the same? Not what I expected


```r
all.equal(c(cov2cor(vv1)),c(fit.us$obj$env$report(fit.us$fit$parfull)$corr[[1]]))
```

```
## [1] TRUE
```

Profiling (experimental/exploratory):


```r
## want $par, not $parfull: do NOT include conditional modes/'b' parameters
ppar <- fit.us$fit$par
length(ppar)
```

```
## [1] 22
```

```r
range(which(names(ppar)=="theta")) ## the last n*(n+1)/2 parameters
```

```
## [1]  2 22
```

```r
## only 1 fixed effect parameter
tt <- tmbprofile(fit.us$obj,2,trace=FALSE)
```


```r
plot(tt)
```

![plot of chunk fit.us.profile.plot](figure/fit.us.profile.plot-1.png)

```r
confint(tt)
```

```
##           lower    upper
## theta 0.3026163 0.390288
```


```r
ppar <- fit.cs$fit$par
length(ppar)
```

```
## [1] 8
```

```r
range(which(names(ppar)=="theta")) ## the last n*(n+1)/2 parameters
```

```
## [1] 2 8
```

```r
## only 1 fixed effect parameter, 1 dispersion parameter
tt2 <- tmbprofile(fit.cs$obj,3,trace=FALSE)
```


```r
plot(tt2)
```

![plot of chunk fit.cs.profile.plot](figure/fit.cs.profile.plot-1.png)
