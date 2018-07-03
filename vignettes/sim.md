---
title: "Simulate from a fitted glmmTMB model"
author: "Mollie Brooks"
date: "2018-06-19"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{simulate}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---
`glmmTMB` has the capability to simulate from a fitted model. These simulations resample random effects from their estimated distribution. In future versions of `glmmTMB`, it may be possible to condition on estimated random effects.




```r
library(glmmTMB)
library(ggplot2); theme_set(theme_bw())
```

Fit a typical model:

```r
data(Owls)
owls_nb1 <- glmmTMB(SiblingNegotiation ~ FoodTreatment*SexParent +
                             (1|Nest)+offset(log(BroodSize)),
                          family = list(family="nbinom1",link="log"),
                          ziformula = ~1, data=Owls)
```

```
## Warning in glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent + (1 | :
## some components missing from 'family': downstream methods may fail
```

```
## Warning in mkTMBStruc(formula, ziformula, dispformula, combForm, mf, fr, :
## specifying 'family' as a plain list is deprecated
```

Then we can simulate from the fitted model with the `simulate.glmmTMB` function. It produces a list of simulated observation vectors, each of which is the same size as the original vector of observations. The default is to only simulate one vector (`nsim=1`) but we still return a list for consistency.


```r
simo=simulate(owls_nb1, seed=1)
Simdat=Owls
Simdat$SiblingNegotiation=simo[[1]]
Simdat=transform(Simdat,  
			NegPerChick = SiblingNegotiation/BroodSize, 
			type="simulated")
Owls$type = "observed"	
Dat=rbind(Owls, Simdat)	
```

Then we can plot the simulated data against the observed data to check if they are similar. 

```r
ggplot(Dat,  aes(NegPerChick, colour=type))+geom_density()+facet_grid(FoodTreatment~SexParent)
```

![plot of chunk plots](figure/plots-1.png)
