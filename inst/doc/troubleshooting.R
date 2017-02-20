## ----echo=FALSE----------------------------------------------------------
library(glmmTMB)

## ----non-pos-def---------------------------------------------------------
zinbm0 = glmmTMB(count~spp + (1|site), zi=~spp, Salamanders, family="nbinom2")

## ----genpois NaN---------------------------------------------------------
m1 = glmmTMB(count~spp + mined + (1|site), zi=~spp + mined, Salamanders, family="genpois")

## ----NA function, eval=FALSE---------------------------------------------
#  Warning in nlminb(start = par, objective = fn, gradient = gr) :
#    NA/NaN function evaluation

## ----Cholmod, eval=FALSE-------------------------------------------------
#  Cholmod warning 'matrix not positive definite'

## ----lgamma, eval=FALSE--------------------------------------------------
#  Warning in f(par, order = order, ...) : value out of range in 'lgamma'

## ----NA gradient, error=TRUE, warning=FALSE------------------------------
dat1 = expand.grid(y=-1:1, rep=1:10)
m1 = glmmTMB(y~1, dat1, family="nbinom2")

## ----gradient length nlminb, eval=FALSE----------------------------------
#  Error in nlminb(start = par, objective = fn, gradient = gr) :
#    gradient function must return a numeric vector of length x

## ----gradient length optimhess, eval=FALSE-------------------------------
#  Error in optimHess(par.fixed, obj$fn, obj$gr):
#    gradient in optim evaluated to length x

