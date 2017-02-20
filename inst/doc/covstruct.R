## ----setup, include=FALSE, message=FALSE---------------------------------
library(knitr)
library(glmmTMB)
library(MASS)
knitr::opts_chunk$set(echo = TRUE, eval=TRUE)
set.seed(1)

## ----echo=TRUE, eval=TRUE------------------------------------------------
n <- 6                                              ## Number of time points
x <- mvrnorm(mu = rep(0,n),
             Sigma = .7 ^ as.matrix(dist(1:n)) )    ## Simulate the process using the MASS package
y <- x + rnorm(n)                                   ## Add measurement noise

## ----echo=TRUE, eval=TRUE------------------------------------------------
times <- factor(1:n)
levels(times)

## ----echo=TRUE, eval=TRUE------------------------------------------------
group <- factor(rep(1,n))

## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  glmmTMB(y ~ ar1(times + 0 | group))

## ----echo=FALSE, eval=TRUE-----------------------------------------------
glmmTMB(y ~ ar1(times + 0 | group))

## ----echo=TRUE, eval=TRUE------------------------------------------------
simGroup <- function(g) {
    x <- mvrnorm(mu = rep(0,n),
             Sigma = .7 ^ as.matrix(dist(1:n)) )    ## Simulate the process
    y <- x + rnorm(n)                               ## Add measurement noise
    times <- factor(1:n)
    group <- factor(rep(g,n))
    data.frame(y, times, group)
}
simGroup(1)

## ----echo=TRUE, eval=TRUE------------------------------------------------
dat <- do.call("rbind", lapply(1:1000, simGroup) )

## ----echo=TRUE, eval=TRUE------------------------------------------------
fit.ar1 <- glmmTMB(y ~ ar1(times + 0 | group), data=dat)
fit.ar1

## ----echo=TRUE, eval=TRUE------------------------------------------------
fit.us <- glmmTMB(y ~ us(times + 0 | group), data=dat, dispformula=~0)
fit.us$sdr$pdHess ## Converged ?

## ----echo=TRUE, eval=TRUE------------------------------------------------
VarCorr(fit.us)

## ----echo=TRUE, eval=TRUE------------------------------------------------
fit.toep <- glmmTMB(y ~ toep(times + 0 | group), data=dat)
fit.toep$sdr$pdHess ## Converged ?

## ----echo=TRUE, eval=TRUE------------------------------------------------
VarCorr(fit.toep)

## ----echo=TRUE, eval=TRUE------------------------------------------------
fit.cs <- glmmTMB(y ~ cs(times + 0 | group), data=dat)
fit.cs$sdr$pdHess ## Converged ?

## ----echo=TRUE, eval=TRUE------------------------------------------------
VarCorr(fit.cs)

## ----echo=TRUE, eval=TRUE------------------------------------------------
anova(fit.ar1, fit.toep, fit.us)

## ----echo=TRUE, eval=TRUE------------------------------------------------
anova(fit.cs, fit.toep)

## ----echo=TRUE, eval=TRUE------------------------------------------------
x <- sample(1:2, 10, replace=TRUE)
y <- sample(1:2, 10, replace=TRUE)

## ----echo=TRUE, eval=TRUE------------------------------------------------
pos <- numFactor(x,y)
pos

## ----echo=TRUE, eval=TRUE------------------------------------------------
parseNumLevels(levels(pos))

## ----echo=TRUE, eval=TRUE------------------------------------------------
dat$times <- numFactor(dat$times)
levels(dat$times)

## ----echo=TRUE, eval=TRUE------------------------------------------------
fit.ou <- glmmTMB(y ~ ou(times + 0 | group), data=dat)
fit.ou$sdr$pdHess ## Converged ?

## ----echo=TRUE, eval=TRUE------------------------------------------------
VarCorr(fit.ou)

## ----echo=TRUE, eval=TRUE------------------------------------------------
fit.mat <- glmmTMB(y ~ mat(times + 0 | group), data=dat, dispformula=~0)
fit.mat$sdr$pdHess ## Converged ?

## ----echo=TRUE, eval=TRUE------------------------------------------------
VarCorr(fit.mat)

## ----echo=TRUE, eval=TRUE------------------------------------------------
fit.gau <- glmmTMB(y ~ gau(times + 0 | group), data=dat, dispformula=~0)
fit.gau$sdr$pdHess ## Converged ?

## ----echo=TRUE, eval=TRUE------------------------------------------------
VarCorr(fit.gau)

## ----echo=TRUE, eval=TRUE------------------------------------------------
fit.exp <- glmmTMB(y ~ exp(times + 0 | group), data=dat)
fit.exp$sdr$pdHess ## Converged ?

## ----echo=TRUE, eval=TRUE------------------------------------------------
VarCorr(fit.exp)

