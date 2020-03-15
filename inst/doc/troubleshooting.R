params <-
list(EVAL = FALSE)

## ----load_lib,echo=FALSE------------------------------------------------------
library(glmmTMB)
knitr::opts_chunk$set(eval = if (isTRUE(exists("params"))) params$EVAL else FALSE)

## ----non-pos-def,cache=TRUE, warning=FALSE------------------------------------
#  zinbm0 = glmmTMB(count~spp + (1|site), zi=~spp, Salamanders, family=nbinom2)

## ----fixef_zinbm0-------------------------------------------------------------
#  fixef(zinbm0)

## ----f_zi2--------------------------------------------------------------------
#  ff <- fixef(zinbm0)$zi
#  round(plogis(c(sppGP=unname(ff[1]),ff[-1]+ff[1])),3)

## ----salfit2,cache=TRUE-------------------------------------------------------
#  Salamanders <- transform(Salamanders, GP=as.numeric(spp=="GP"))
#  zinbm0_A = update(zinbm0, ziformula=~GP)

## ----salfit2_coef,cache=TRUE--------------------------------------------------
#  fixef(zinbm0_A)[["zi"]]

## ----salfit3,cache=TRUE-------------------------------------------------------
#  zinbm0_B = update(zinbm0, ziformula=~(1|spp))
#  fixef(zinbm0_B)[["zi"]]
#  VarCorr(zinbm0_B)

## ----zinbm1,cache=TRUE--------------------------------------------------------
#  zinbm1 = glmmTMB(count~spp + (1|site), zi=~mined, Salamanders, family=nbinom2)
#  fixef(zinbm1)[["zi"]]

## ----zinbm1_confint,cache=TRUE------------------------------------------------
#  ## at present we need to specify the parameter by number; for
#  ##  extreme cases need to specify the parameter range
#  ## (not sure why the upper bound needs to be so high ... ?)
#  cc = confint(zinbm1,method="uniroot",parm=9, parm.range=c(-20,20))
#  print(cc)

## ----fatfiberglmm-------------------------------------------------------------
#  ## data taken from gamlss.data:plasma, originally
#  ## http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/plasma.html
#  load(system.file("vignette_data","plasma.rda", package="glmmTMB"))
#  m4.1 <- glm(calories ~ fat*fiber, family = Gamma(link = "log"), data = plasma)
#  m4.2 <- glmmTMB(calories ~ fat*fiber, family = Gamma(link = "log"), data = plasma)
#  ps  <- transform(plasma,fat=scale(fat,center=FALSE),fiber=scale(fiber,center=FALSE))
#  m4.3 <- update(m4.2, data=ps)
#  ## scaling factor for back-transforming standard deviations
#  ss <- c(1,
#          fatsc <- 1/attr(ps$fat,"scaled:scale"),
#          fibsc <- 1/attr(ps$fiber,"scaled:scale"),
#          fatsc*fibsc)
#  ## combine SEs, suppressing the warning from the unscaled model
#  s_vals <- cbind(glm=sqrt(diag(vcov(m4.1))),
#                  glmmTMB_unsc=suppressWarnings(sqrt(diag(vcov(m4.2)$cond))),
#                  glmmTMB_sc=sqrt(diag(vcov(m4.3)$cond))*ss)
#  print(s_vals,digits=3)

## ----load_ss_ex---------------------------------------------------------------
#  load(system.file("vignette_data","troubleshooting.rda",package="glmmTMB"))

## ----ss_ex_mod1---------------------------------------------------------------
#  summary(mod1)

## ----diagnose_vcov------------------------------------------------------------
#  diagnose_vcov <- function(model, tol=1e-5, digits=2, analyze_hessian=FALSE) {
#      vv <- vcov(model, full=TRUE)
#      nn <- rownames(vv)
#      if (!all(is.finite(vv))) {
#          if (missing(analyze_hessian)) warning("analyzing Hessian, not vcov")
#          if (!analyze_hessian) stop("can't analyze vcov")
#          analyze_hessian <- TRUE
#      }
#      if (analyze_hessian) {
#          par.fixed <- model$obj$env$last.par.best
#          r <- model$obj$env$random
#          if (!is.null(r)) par.fixed <- par.fixed[-r]
#          vv <- optimHess(par.fixed, fn=model$obj$fn, gr=model$obj$gr)
#          ## note vv is now HESSIAN, not vcov
#      }
#      ee <- eigen(vv)
#      if (all(ee$values>tol)) {message("var-cov matrix OK"); return(invisible(NULL))}
#      ## find negative or small-positive eigenvalues (flat/wrong curvature)
#      bad_evals <- which(ee$values<tol)
#      ## order worst to best
#      bad_evals <- bad_evals[order(-ee$values[bad_evals])]
#      ret <- lapply(bad_evals,
#                    function(i) {
#                        ## extract loadings
#                        v <- setNames(ee$vectors[,i], nn)
#                        ## order in decreasing magnitude & round
#                        list(val=ee$values[i],vec=round(v[order(-abs(v))],digits))
#                    })
#      return(ret)
#  }

## ----diag_1-------------------------------------------------------------------
#  (d1 <- diagnose_vcov(mod1))

## ----ss_mod2_up, eval=FALSE---------------------------------------------------
#  mod2 <- update(mod1, ziformula=~0)

## ----ss_mod2------------------------------------------------------------------
#  summary(mod2)

## ----ss_diag2-----------------------------------------------------------------
#  diagnose_vcov(mod2)

## ----mod3_up, eval=FALSE------------------------------------------------------
#  mod3 <- update(mod2, family=poisson)

## ----ss_mod3------------------------------------------------------------------
#  summary(mod3)

## ----ss_diag3-----------------------------------------------------------------
#  diagnose_vcov(mod3)

## ----checkhess----------------------------------------------------------------
#  mod3$sdr$pdHess					

## ----genpois_NaN,cache=TRUE---------------------------------------------------
#  m1 = glmmTMB(count~spp + mined + (1|site), zi=~spp + mined, Salamanders, family=genpois)

## ----NA gradient, error=TRUE, warning=FALSE-----------------------------------
#  dat1 = expand.grid(y=-1:1, rep=1:10)
#  m1 = glmmTMB(y~1, dat1, family=nbinom2)

