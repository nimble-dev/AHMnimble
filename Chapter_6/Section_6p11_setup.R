#---(begin AHMnimble header)---
# This file contains code adapted from the file R_BUGS_code_AHM_Vol_1_20170519.R
# available at https://www.mbr-pwrc.usgs.gov/pubanalysis/keryroylebook/, with
# permission of the authors.  It may have been modified so that examples originally
# written to work with JAGS or WinBUGS will work with NIMBLE.  More information
# about NIMBLE can be found at https://r-nimble.org, https://github.com/nimble-dev/nimble,
# and https://cran.r-project.org/web/packages/nimble/index.html.
#
# The file  R_BUGS_code_AHM_Vol_1_20170519.R contains the following header:
# -----(begin R_BUGS_code_AHM_Vol_1_20170519.R header)-----
# =========================================================================
#
#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#
#   Marc Kéry & J. Andy Royle
#
#   *** This is the text file with all R and BUGS code from the book ***
#
#   Created 2 Dec 2015 based on draft from 21 Oct 2015
#
# =========================================================================
### Last change: 19 May 2017 by Mike Meredith
# Incorporated errata up to 19 May 2017
# Code updated to implement new names and changes to functions in AHMbook package 0.1.4
# Use built-in data sets instead of .csv files
# In Chapter 11, replaced 'Y', 'Ysum' and 'Yaug' with lower case 'y', 'ysum' and 'yaug'
#  to match the code in the printed book.
#
# -----(end R_BUGS_code_AHM_Vol_1_20170519.R header)-----

# This file was created by:
#
# Jacob Levine and Perry de Valpine
#
#---(end AHMnimble header)---
source("Chapter_6_setup_nimble.R")
# 6.11.1 Bayesian fitting of the basic ZIP N-mixture model
# ------------------------------------------------------------------------
# From section 6.9.2 - for data creation/organization:

library(AHMbook)
## Code modified to use the SwissTits data set included in the AHMbook package
data(SwissTits)
str(SwissTits)
SwissTits$species  # Available species

# Select Great tit and covariate data from 2013 and
#   drop 4 sites not surveyed in 2013
y0 <- SwissTits$counts[, , '2013', 'Great tit']
( NA.sites <- which(rowSums(is.na(y0)) == 3) ) # Unsurveyed sites
y <- y0[-NA.sites, ]                 # Drop them from the count data
tits <- SwissTits$sites[-NA.sites, ] # Also drop from the site covariates
str(y)  # Check the matrix of count data
# Get date and duration data for 2013, without the NA.sites rows:
date <- SwissTits$date[-NA.sites, , '2013']
dur <- SwissTits$dur[-NA.sites, , '2013']

# Plot observed data: counts vs survey date (Fig. 6-9)
matplot(t(date), t(y), type = "l", lwd = 3, lty = 1, frame = F, xlab = "Survey data (1 = April 1)", ylab = "Count of Great Tits")

# Load unmarked, create unmarked data frame and inspect result
library(unmarked)
time <- matrix(rep(as.character(1:3), nrow(y)), ncol = 3, byrow = TRUE)
umf <- unmarkedFramePCount(y = y,
                           siteCovs=data.frame(elev=scale(tits[,"elev"]), forest=scale(tits[,"forest"]), iLength=1/tits[,"rlength"]),
                           obsCovs=list(time = time, date = scale(date), dur = scale(dur)))
summary(umf)                            # Summarize unmarked data frame
summary(apply(y, 1, max, na.rm = TRUE)) # Summarize max counts

elev <- umf@siteCovs$elev   ;   elev2 <- elev^2
forest <- umf@siteCovs$forest   ;   forest2 <- forest^2
date <- matrix(umf@obsCovs$date, ncol = 3, byrow = TRUE)
dur <- matrix(umf@obsCovs$dur, ncol = 3, byrow = TRUE)
date[is.na(date)] <- 0   ;   date2 <- date^2
dur[is.na(dur)] <- 0   ;   dur2 <- dur^2
iRoute <- umf@siteCovs$iLength

# Design matrix for abundance model (no intercept)
lamDM <- model.matrix(~ elev + elev2 + forest + forest2 + elev:forest + elev:forest2 + iRoute)[,-1]

require(nimble)
# Specify model in BUGS language
# "ZIPNmix.txt"
Section6p11_code <- nimbleCode( {
    
    # Specify priors
    # zero-inflation/suitability
    phi ~ dunif(0,1)          # proportion of suitable sites (probability of being not a structural 0)
    theta <- 1-phi            # zero-inflation (proportion of unsuitable)
    ltheta <- logit(theta)
    
    # abundance
    beta0 ~ dnorm(0, 0.1)     # log(lambda) intercept
    for(k in 1:7){            # Regression params in lambda
      beta[k] ~ dnorm(0, 1)
    }
    tau.lam <- pow(sd.lam, -2)
    sd.lam ~ dunif(0, 2)      # site heterogeneity in lambda
    
    # detection
    for(j in 1:3){
      alpha0[j] <- logit(mean.p[j])
      mean.p[j] ~ dunif(0, 1) # p intercept for occasions 1-3
      }
    for(k in 1:13){           # Regression params in p
      alpha[k] ~ dnorm(0, 1)
      }
    tau.p.site <- pow(sd.p.site, -2)
    sd.p.site ~ dunif(0, 2)   # site heterogeneity in p
    tau.p.survey <- pow(sd.p.survey, -2)
    sd.p.survey ~ dunif(0, 2) # site-survey heterogeneity in p
    
    # ZIP model for abundance
    for (i in 1:nsite){
      a[i] ~ dbern(phi)
      eps.lam[i] ~ dnorm(0, tau.lam)       # Random site effects in log(abundance)
      loglam[i] <- beta0 + inprod(beta[1:7], lamDM[i, 1:7]) + eps.lam[i] * hlam.on
      loglam.lim[i] <- min(250, max(-250, loglam[i]))  # �Stabilize� log
      lam[i] <- exp(loglam.lim[i])
      mu.poisson[i] <- a[i] * lam[i]
      N[i] ~ dpois(mu.poisson[i])
    }
    
    # Measurement error model
    for (i in 1:nsite){
      eps.p.site[i] ~ dnorm(0, tau.p.site) # Random site effects in logit(p)
      for (j in 1:nrep){
        y[i,j] ~ dbin(p[i,j], N[i])
        p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
        lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  # �Stabilize� logit
        lp[i,j] <- alpha0[j] + alpha[1] * elev[i] + alpha[2] * elev2[i] +
          alpha[3] * date[i,j] + alpha[4] * date2[i,j] +
          alpha[5] * dur[i,j] + alpha[6] * dur2[i,j] +
          alpha[7] * elev[i] * date[i,j] + alpha[8] * elev2[i] * date[i,j] +
          alpha[9] * elev[i] * dur[i,j] + alpha[10] * elev[i] * dur2[i,j] +
          alpha[11] * elev2[i] * dur[i,j] + alpha[12] * date[i,j] * dur[i,j] +
          alpha[13] * date[i,j] * dur2[i,j] +
          eps.p.site[i] * hp.site.on + eps.p.survey[i,j] * hp.survey.on
          eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) # Random site-survey effects
    }
    }
    # Posterior predictive distributions of chi2 discrepancy
    for (i in 1:nsite) {
      for (j in 1:nrep) {
        y.sim[i,j] ~ dbin(p[i,j], N[i]) # Create new data set under model
        e.count[i,j] <- N[i] * p[i,j]   # Expected datum
        # Chi-square discrepancy for the actual data
        chi2.actual[i,j] <- pow((y[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
        # Chi-square discrepancy for the simulated ('perfect') data
        chi2.sim[i,j] <- pow((y.sim[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
    # Add small value e to denominator to avoid division by zero
    }
    }
    # Add up individual chi2 values for overall fit statistic
    fit.actual <- sum(chi2.actual[1:263, 1:3])  # Fit statistic for actual data set
    fit.sim <- sum(chi2.sim[1:263, 1:3])        # Fit statistic for a fitting model
    bpv <- step(fit.sim-fit.actual)    # Bayesian p-value
    c.hat <- fit.actual/fit.sim        # c-hat estimate
    
    # Derived parameters: Total abundance at 263 sampled sites
    Ntotal263 <- sum(N[1:263])
    }
)

# Initial values
Nst <- apply(y, 1, max, na.rm = T) + 1
Nst[is.na(Nst)] <- round(mean(y, na.rm = TRUE))
Nst[Nst == "-Inf"] <- round(mean(y, na.rm = TRUE))
inits <- function(){ list(N = Nst, 
                          beta0 = 0, 
                          mean.p = rep(0.5, 3), 
                          beta = runif(7, 0, 0), 
                          alpha = runif(13, 0, 0)
                          )}

# Bundle data and choose to fit simple ZIP model (model 1)
win.data1 <- list(y = y, 
                  lamDM = lamDM, 
                  elev = elev, 
                  date = date, 
                  dur = dur, 
                  elev2 = elev2,
                  date2 = date2, 
                  dur2 = dur2, 
                  e = 1e-06, 
                  hlam.on = 0, 
                  hp.site.on = 0,
                  hp.survey.on = 0,
                  nsite = 263,
                  nrep = 3)

SGT_data1 <- win.data1

## added for nimble:
SGT_inits_full<- function(){
    ans <- list(N = Nst, 
                beta0 = 0, 
                mean.p = rep(0.5, 3), 
                beta = runif(7, 0, 0), 
                alpha = runif(13, 0, 0)
              , phi = 0.5
              , sd.lam = 1
              , sd.p.site = 1
              , sd.p.survey = 1
              , a = rbinom(SGT_data1$nsite, prob = 0.5, size = 1)
              , eps.lam = rnorm(SGT_data1$nsite)
              , eps.p.site = rnorm(SGT_data1$nsite)
              , eps.p.survey = matrix(rnorm(SGT_data1$nsite * SGT_data1$nrep),
                                      nrow = SGT_data1$nsite)
                )
    ans$a[ ans$N > 0 ] <- 1
    ans
}
