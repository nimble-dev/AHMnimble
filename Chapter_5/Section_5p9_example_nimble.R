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
source("Chapter_5_setup.R")
# 5.9 Poisson generalized linear model (Poisson GLM)
# --------------------------------------------------

# Bundle data
win.data <- list(Cmax = Cmax, M = length(Cmax), elev = elev, facFor = facFor, e = 0.0001)

require(nimble)
# Specify model in BUGS language
# Code for file "Poisson_GLM.txt" from AHM code.
section5p9Code <- nimbleCode( {
    
    # Priors
    for(k in 1:4){
    alpha[k] ~ dnorm(0, 1.0E-06)       # Prior for intercepts
    beta[k] ~ dnorm(0, 1.0E-06)        # Prior for slopes
    }
    
    # Likelihood
    for (i in 1:M){
    Cmax[i] ~ dpois(lambda[i])         # note no variance parameter
    log(lambda[i]) <- alpha[facFor[i]] + beta[facFor[i]] * elev[i]
    resi[i] <- (Cmax[i]-lambda[i]) / (sqrt(lambda[i])+e)   # Pearson resi
    }
  }
)


# Initial values
inits <- function() list(alpha = rnorm(4,,3), beta = rnorm(4,,3))

# Parameters monitored
params <- c("alpha", "beta", "lambda", "resi")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Build, configure and compile model and MCMC
out5 <- nimbleMCMC(code = section5p9Code,
                   constants = win.data,
                   inits = inits,
                   monitors = params,
                   nburnin = nb,
                   niter = ni)

require(coda)
out5.mcmc <- as.mcmc(out5)
require(mcmcplots)
mcmcplot(out5.mcmc)

## Other AHM book plots skipped

