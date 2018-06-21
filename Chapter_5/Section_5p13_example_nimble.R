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
# 5.13 Random-effects Poisson GLM (Poisson GLMM)
# ----------------------------------------------

# Bundle data
win.data <- list(C = C, M = nrow(C), J = ncol(C), elev = elev, forest = forest, elev.forest = elev * forest, wind = wind)
require(nimble)
# Specify model in BUGS language
# Code for file "RE.Poisson.txt" from AHM code
section5p13_code <- nimbleCode( {
    
    # Priors
    mu.alpha ~ dnorm(0, 0.001)                # Mean hyperparam
    tau.alpha <- pow(sd.alpha, -2)
    sd.alpha ~ dunif(0, 10)                   # sd hyperparam
    for(k in 1:4){
    alpha[k] ~ dunif(-10, 10)              # Regression params
    }
    
    # Likelihood
    for (i in 1:M){
    alpha0[i] ~ dnorm(mu.alpha, tau.alpha) # Random effects and hyperparams
    re0[i] <- alpha0[i] - mu.alpha         # zero-centered random effects
    for(j in 1:J){
    C[i,j] ~ dpois(lambda[i,j])
    log(lambda[i,j]) <- alpha0[i] + alpha[1] * elev[i] + alpha[2] * forest[i] + alpha[3] * elev.forest[i] + alpha[4] * wind[i,j]
    }
  }
}
)
# Other model run preparations
inits <- function() list(alpha0 = rnorm(M), alpha = rnorm(4)) # Inits
params <- c("mu.alpha", "sd.alpha", "alpha0", "alpha", "re0") # Params
ni <- 30000 ; nt <- 25 ; nb <- 5000 ; nc <- 3                 # MCMC settings

# Call nimble from R:
out8 <- nimbleMCMC(code = section5p13_code,
                   constants = win.data,
                   inits = inits,
                   monitors = params,
                   nburnin = nb,
                   niter = ni)
library(coda)
out8.mcmc <- as.mcmc(out8)
dim(out8.mcmc)
## This is quite large, so no figures are produced.
