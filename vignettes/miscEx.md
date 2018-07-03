---
title: "Miscellaneous examples"
date: "2018-06-19"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{miscellaneous examples}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



## Beta dispersion model


```r
set.seed(1001)
N <- 1000
mean_pars <- c(1,2)
disp_pars <- c(1,2)
dd <- data.frame(x=rnorm(N))
m <- plogis(mean_pars[1]+mean_pars[2]*dd$x)
d <- exp(disp_pars[1]+disp_pars[2]*dd$x)
dd$y <- rbeta(N,shape1=m*d,shape2=(1-m)*d)
```

Fit models:


```r
## location only
m1 <- glmmTMB(y~x,
              family=beta_family(),
              data=dd)
## add model for dispersion
m2 <- update(m1,dispformula=~x)
```

Fixed effects look close to theoretical values:

```r
fixef(m2)
```

```
## 
## Conditional model:
## (Intercept)            x  
##       1.005        2.013  
## 
## Dispersion model:
## (Intercept)            x  
##       1.064        1.962
```

AIC is insanely much better for the model with dispersion varying:

```r
bbmle::AICtab(m1,m2)
```

```
##    dAIC   df
## m2    0.0 4 
## m1 1491.6 3
```
