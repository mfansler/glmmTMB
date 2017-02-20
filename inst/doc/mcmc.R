## ----setup, include=FALSE, message=FALSE---------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE)

## ----libs,message=FALSE--------------------------------------------------
library(glmmTMB)
library(MCMCpack)
library(coda)
library(lattice)
library(ggplot2); theme_set(theme_bw())
library(reshape2)

## ----fit1----------------------------------------------------------------
data("sleepstudy",package="lme4")
fm1 <- glmmTMB(Reaction ~ Days + (Days|Subject),
               sleepstudy)

## ----mcmcsetup-----------------------------------------------------------
## FIXME: better way for user to extract full coefs?
rawcoef <- with(fm1$obj$env,last.par[-random])
stopifnot(all.equal(c(fm1$obj$fn(rawcoef)),
                    -c(logLik(fm1)),
          tolerance=1e-7))
## log-likelihood function 
## (MCMCmetrop1R wants *positive* log-lik)
fn <- function(x) -fm1$obj$fn(x)
V <- vcov(fm1,full=TRUE)

## ----do_mcmc,echo=FALSE--------------------------------------------------
fn <- system.file("vignette_data","mcmc.rda",package="glmmTMB")
if (length(fn)==0) {
    s1 <- system.time(m1 <- MCMCmetrop1R(fn,rawcoef,V=V))
    save("s1","m1",file="mcmc.rda")
} else {
    load(fn)
}

## ----show_mcmc,eval=FALSE------------------------------------------------
#  s1 <- system.time(m1 <- MCMCmetrop1R(fn,rawcoef,V=V))

## ------------------------------------------------------------------------
colnames(m1) <- c(names(fixef(fm1)[[1]]),
                        "log(sigma)",
                  c("log(sd_Intercept)","log(sd_Days)","cor"))

## ----traceplot,fig.width=7-----------------------------------------------
xyplot(m1,layout=c(2,3))

## ----effsize-------------------------------------------------------------
print(effectiveSize(m1),digits=3)

## ----violins,echo=FALSE--------------------------------------------------
ggplot(melt(as.matrix(m1[,-1])),aes(x=Var2,y=value))+
         geom_violin(fill="gray")+coord_flip()+labs(x="")

## ----echo=FALSE,eval=FALSE-----------------------------------------------
#  s2 <- system.time(m2 <- MCMCmetrop1R(fn,rawcoef,V=V,
#                                       tune=c(1,1,1,1,1,100)))

