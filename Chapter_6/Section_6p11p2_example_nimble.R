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
source("Section_6p11_setup.R")
# 6.11.2.1 Accounting for overdispersion at multiple scales
# ------------------------------------------------------------------------
# MCMC settings
ni <- 10^6    ;    nt <- 80   ;    nb <- 200000    ;    nc <- 3
ni <- 100; nb <- 10; nt <- 2;
# This section runs MCMC for a variety of models established by different
# values (0 or 1) of hlam.on, hp.site.on, and hb.survey.on.
# In NIMBLE, the values of these can be changed in the same model object,
# and then the same MCMC can be re-run.  To do so, we need to more carefully
# separate data from constants.  Otherwise hlam.on, hp.site.on, and hp.survey.on
# will be treated as constants and will not be modifiable.
Section6p11_constants <- win.data1[c('e', 'nsite', 'nrep')]
Section6p11_data <- win.data1
Section6p11_data[c('e','nsite','nrep')] <- NULL

Section6p11_model <- nimbleModel(Section6p11_code, data = Section6p11_data, constants = Section6p11_constants, inits = SGT_inits_full())
Section6p11_MCMC <- buildMCMC(Section6p11_model, thin = nt)
cSection6p11_model <- compileNimble(Section6p11_model)
cSection6p11_MCMC <- compileNimble(Section6p11_MCMC, project = Section6p11_model)

# Bundle data and select model 2
## win.data2 <- list(y = y, nsite = nrow(y), nrep = ncol(y), 
##                   lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2, 
##                   date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 1, hp.site.on = 0,
##                   hp.survey.on = 0)
cSection6p11_model$hlam.on <- 1
cSection6p11_model$hp.site.on <- 0
cSection6p11_model$hp.survey.on <- 0
out2 <- runMCMC(cSection6p11_MCMC, niter = ni, nburnin = nb, nchains = nc)
require(mcmcplots)
mcmcplot(out2)


# Bundle data and select model 3
## win.data3 <- list(y = y, nsite = nrow(y), nrep = ncol(y), 
##                   lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2, 
##                   date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 0, hp.site.on = 1,
##                   hp.survey.on = 0)

cSection6p11_model$hlam.on <- 0
cSection6p11_model$hp.site.on <- 1
cSection6p11_model$hp.survey.on <- 0
out3 <- runMCMC(cSection6p11_MCMC, niter = ni, nburnin = nb, nchains = nc)
require(mcmcplots)
mcmcplot(out3)


# Bundle data and select model 4
## win.data4 <- list(y = y, nsite = nrow(y), nrep = ncol(y), 
##                   lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2, 
##                   date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 0, hp.site.on = 0,
##                   hp.survey.on = 1)

cSection6p11_model$hlam.on <- 0
cSection6p11_model$hp.site.on <- 0
cSection6p11_model$hp.survey.on <- 1
out4 <- runMCMC(cSection6p11_MCMC, niter = ni, nburnin = nb, nchains = nc)
require(mcmcplots)
mcmcplot(out4)


# Bundle data and select model 5
## win.data5 <- list(y = y, nsite = nrow(y), nrep = ncol(y), 
##                   lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2, 
##                   date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 1, hp.site.on = 1,
##                   hp.survey.on = 0)
cSection6p11_model$hlam.on <- 1
cSection6p11_model$hp.site.on <- 1
cSection6p11_model$hp.survey.on <- 0
out5 <- runMCMC(cSection6p11_MCMC, niter = ni, nburnin = nb, nchains = nc)
require(mcmcplots)
mcmcplot(out5)


# Bundle data and select model 6
## win.data6 <- list(y = y, nsite = nrow(y), nrep = ncol(y), 
##                   lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2, 
##                   date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 1, hp.site.on = 0,
##                   hp.survey.on = 1)
cSection6p11_model$hlam.on <- 1
cSection6p11_model$hp.site.on <- 0
cSection6p11_model$hp.survey.on <- 1
out6 <- runMCMC(cSection6p11_MCMC, niter = ni, nburnin = nb, nchains = nc)
require(mcmcplots)
mcmcplot(out6)


# Bundle data and select model 7
## win.data7 <- list(y = y, nsite = nrow(y), nrep = ncol(y), 
##                   lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2, 
##                   date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 1, hp.site.on = 1,
##                   hp.survey.on = 1)
cSection6p11_model$hlam.on <- 1
cSection6p11_model$hp.site.on <- 1
cSection6p11_model$hp.survey.on <- 1
out7 <- runMCMC(cSection6p11_MCMC, niter = ni, nburnin = nb, nchains = nc)
require(mcmcplots)
mcmcplot(out7)



# 6.11.2.2 Linear modeling of a variance in the N-mixture model
# ------------------------------------------------------------------------
# Bundle and summarize data set
str(win.data8 <- list(y = y, nsite = nrow(y), nrep = ncol(y), 
                       lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2, 
                       date2 = date2, dur2 = dur2) )

# Specify model in BUGS language
# This corresponds to "Nmix.special.txt" in original book code
section6p11p2_code <- nimbleCode( {
    
    # Specify priors
    # abundance
      beta0 ~ dnorm(0, 0.1)       # log(lambda) intercept
        for(k in 1:7){            # Regression params in lambda
          beta[k] ~ dnorm(0, 1)
        }
    # Model for unexplained variance in lambda among sites
      for (i in 1:nsite){
        tau.lam[i] <- 1/var.lam[i]
        log(var.lam[i]) <- alpha.var.lam + beta.var.lam * elev[i]
      }
    # Priors for intercept and slope of linear model for variance
      alpha.var.lam ~ dunif(-1, 1)
      beta.var.lam ~ dunif(0, 3)
      
    # detection
      for(j in 1:3){
        alpha0[j] <- logit(mean.p[j])
        mean.p[j] ~ dunif(0, 1)# p intercept for occasions 1-3
      }
      for(k in 1:13){           # Regression params in p
        alpha[k] ~ dnorm(0, 1)
      }
      tau.p.survey <- pow(sd.p.survey, -2)
      sd.p.survey ~ dunif(0, 1) # site-survey heterogeneity in p
    
    # Poisson-lognormal model for abundance
      for (i in 1:nsite){
        eps.lam[i] ~ dnorm(0, tau.lam[i]) # Random site effects in log(abundance)
        loglam[i] <- beta0 + inprod(beta[1:7], lamDM[i,]) + eps.lam[i]
        loglam.lim[i] <- min(250, max(-250, loglam[i]))  # �Stabilize� log
        mu.poisson[i] <- exp(loglam.lim[i])
        N[i] ~ dpois(mu.poisson[i])
      }
    
    # Binomial measurement error model with extra-binomial dispersion
    for (i in 1:nsite){
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
        alpha[13] * date[i,j] * dur2[i,j] + eps.p.survey[i,j]
        eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) # Random site-survey effects
      }
    }
  }
 )


# Initial values
Nst <- apply(y, 1, max, na.rm = T) + 1
Nst[is.na(Nst)] <- round(mean(y, na.rm = TRUE))
Nst[Nst == "-Inf"] <- round(mean(y, na.rm = TRUE))
inits <- function(){ list(N = Nst, 
                          beta0 = 0, 
                          mean.p = rep(0.5,3), 
                          beta = runif(7, 0,0), 
                          alpha = runif(13, 0,0), 
                          alpha.var.lam = 0, 
                          beta.var.lam = 1.5, 
                          sd.p.survey = 0.3)}

# Parameters monitored
params <- c("beta0", "beta", "alpha.var.lam", "beta.var.lam", "alpha0", "mean.p", "alpha", "sd.p.survey")

# MCMC settings
ni <- 180000    ;    nt <- 100    ;    nb <- 10000    ;    nc <- 3

# Call nimble from R (ART 374 min) and summarize posteriors
out8 <- nimbleMCMC(code = section6p11p2_code, 
                   constants = win.data8, 
                   inits = SGT_inits_full,
                   monitors = params,
                   nburnin = nb, 
                   niter = ni,
                   thin = nt,
                   samplesAsCodaMCMC = TRUE)

library(coda)
colnames(out8)
library(mcmcplots)
mcmcplot(out8)

# Compare MCMC performance across packages:
params_of_interest = c("beta0", "beta", "alpha.var.lam", "beta.var.lam", "alpha0", "mean.p", "alpha", "sd.p.survey")
Section6p11p2_compare <- compareMCMCs(
  modelInfo = list(
    code = section6p11p2_code,
    data = win.data8,
    inits = SGT_inits_full()),
  monitors = params_of_interest,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = nb,
  thin = nt,
  niter = ni
)

make_MCMC_comparison_pages(Section6p11p2_compare, modelNames = "Section6p11p2", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section6p11p2.html"))

