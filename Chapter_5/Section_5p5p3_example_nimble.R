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
source("Chapter_5_setup")

# 5.5.3 Missing values in a covariate
# ------------------------------------------------------------------------
## Shoot 'holes' in the covariate data
ele <- elev          # copy of elevation covariate
ele[1:10] <- NA      # missing values in covariate elevation

## Bundle data: feed new 'ele' into 'elev' covariate inside of BUGS model
win.data <- list(Cmean = Cmean, M = length(Cmean), elev = ele, forest = forest)

## Specify model in NIMBLE language (method 2 - mean imputation):

require(nimble)

# Code for file "missing_cov_imputation_model_1.txt"
# from original AHM code.
Section5p5p3_code <- nimbleCode({
    
    # Priors
    alpha0 ~ dnorm(0, 1.0E-06)           # Prior for intercept
    alpha1 ~ dnorm(0, 1.0E-06)           # Prior for slope of elev
    alpha2 ~ dnorm(0, 1.0E-06)           # Prior for slope of forest
    alpha3 ~ dnorm(0, 1.0E-06)           # Prior for slope of interaction
    tau <- pow(sd, -2)
    sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale
    
    # Likelihood
    for (i in 1:M){
    Cmean[i] ~ dnorm(mu[i], tau)      # precision tau = 1 / variance
    mu[i] <- alpha0 + alpha1 * elev[i] + alpha2 * forest[i] + alpha3 * elev[i] * forest[i]
    }
    
    # Model for missing covariates
    for (i in 1:M){
    elev[i] ~ dnorm(0, 0.001)
    }
})

## Initial values
inits <- function() list(alpha0 = rnorm(1,,10), alpha1 = rnorm(1,,10), alpha2 = rnorm(1,,10), alpha3 = rnorm(1,,10))

## Parameters monitored
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "elev")

## MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

out1.3 <- nimbleMCMC(code = Section5p5p3_code,
                     constants = win.data,
                     inits = inits,
                     monitors = params,
                     nburnin = nb,
                     niter = ni)

require(coda)
out1.3.mcmc <- as.mcmc(out1.3)
require(mcmcplots)
mcmcplot(out1.3.mcmc[,1:14])

## Specify model in BUGS language (method 3)
## Code for file "missing_cov_imputation_model_2.txt" in original book code.
section5p5p3_code.3 <-  nimbleCode({
    
    # Priors
    alpha0 ~ dnorm(0, 1.0E-06)           # Prior for intercept
    alpha1 ~ dnorm(0, 1.0E-06)           # Prior for slope of elev
    alpha2 ~ dnorm(0, 1.0E-06)           # Prior for slope of forest
    alpha3 ~ dnorm(0, 1.0E-06)           # Prior for slope of interaction
    tau <- pow(sd, -2)
    sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale
    
    # Likelihood
    for (i in 1:M){
    Cmean[i] ~ dnorm(mu[i], tau)      # precision tau = 1 / variance
    mu[i] <- alpha0 + alpha1 * elev[i] + alpha2 * forest[i] + alpha3 * elev[i] * forest[i]
    }
    
    # Covariate mean as a model for missing covariates
    for (i in 1:M){
    elev[i] ~ dnorm(mu.elev, tau.elev)    # Assume elevation normally distributed
    }
    mu.elev ~ dnorm(0, 0.0001)
    tau.elev <- pow(sd.elev, -2)
    sd.elev ~ dunif(0, 100)
    })

## Initial values
inits <- function() list(alpha0 = rnorm(1,,10), alpha1 = rnorm(1,,10), alpha2 = rnorm(1,,10), alpha3 = rnorm(1,,10))

# # Parameters monitored
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "elev", "mu.elev", "sd.elev")

## MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

## Build, Configure and Compile MCMC:

out1.4 <- nimbleMCMC(code = section5p5p3_code.3,
                     constants = win.data,
                     inits = inits,
                     monitors = params,
                     nburnin = nb,
                     niter = ni)

out1.4.mcmc <- as.mcmc(out1.4)
mcmcplot(out1.4.mcmc[,1:14])
