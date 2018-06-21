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
# 6.11.1 Bayesian fitting of the basic ZIP N-mixture model
# ------------------------------------------------------------------------

# Note that this model contains "toggle" parameters that switch parts of the
# model on and off.  In this case, the three random effects
# are all turned off (hlam.on = hp.survey.on = hp.site.on = 0).
# This means that the random effects and their standard deviations 
# (sd.lam, sd.p.site, sd.p.survey) will just follow their priors.
# Nevertheless, they are being sampled, and this results in a substantial
# amount of wasted computation time.
# In the comparisons later, we will turn off the relevant samplers in NIMBLE.

# Parameters monitored
params <- c("theta", "ltheta", "phi", "beta0", "beta", "alpha0", "mean.p", "alpha", "fit.actual", "fit.sim", "bpv", "c.hat", "Ntotal263")

# MCMC settings
ni <- 50000    ;    nt <- 4    ;    nb <- 10000    ;    nc <- 3

# Call nimble from R and summarize posteriors:
out1 <- nimbleMCMC(code = Section6p11_code, 
                   constants = win.data1, 
                   inits = SGT_inits_full, 
                   monitors = params,
                   nburnin = nb, 
                   niter = ni,
                   thin = nt,
                   samplesAsCodaMCMC = TRUE)

## Note that we provided more complete intitial values for NIMBLE than was
## provied in the original book code for JAGS.

require(coda)
require(mcmcplots)
colnames(out1)
## NIMBLE: We choose to monitor the actual parameters, not derived quantities, for comparison purposes.
params_of_interest <- c('Ntotal263', 'alpha[1]', 'alpha[6]', 'alpha0[1]', 'alpha0[3]', 'beta[1]', 'beta[6]', 'beta0','ltheta', 'mean.p[1]', 'phi', 'theta')
mcmcplot(out1[, params_of_interest])

set.seed(1)
initsValues = SGT_inits_full()

## For this example, we will set up some different sampler configurations in NIMBLE:
## 1. A configuration that does not waste time sampling nodes that are turned "off" by toggle variables.
## 2. A configuration like (1) and also samples some nodes by Automated Factor Slice Samplers.
## 3. A configuration like (1) and also samples some nodes by Adaptive Random-Walk Metropolis-Hastings block samplers.

removeUnusedSamplers <- function(MCMCconf) {
    MCMCconf$removeSamplers("eps.lam")
    MCMCconf$removeSamplers("eps.p.site")
    MCMCconf$removeSamplers("eps.p.survey")
    MCMCconf$removeSamplers("sd.lam")
    MCMCconf$removeSamplers("sd.p.site")
    MCMCconf$removeSamplers("sd.p.survey")
    MCMCconf
}

assignAFSS <- function(MCMCconf) {
    MCMCconf$removeSamplers(c("beta0","beta"))
    MCMCconf$addSampler(c("beta0","beta"), type = "AF_slice")
    MCMCconf$removeSamplers(c("mean.p", "alpha"))
    MCMCconf$addSampler(c("mean.p", "alpha"), type = "AF_slice")
    MCMCconf
}

assignRWB <- function(MCMCconf) {
    MCMCconf$removeSamplers(c("beta0","beta"))
    MCMCconf$addSampler(c("beta0","beta"), type = "RW_block")
    MCMCconf$removeSamplers(c("mean.p", "alpha"))
    MCMCconf$addSampler(c("mean.p", "alpha"), type = "RW_block")
    MCMCconf
}

MCMCdefs = list(
    nimble2 = quote({
        removeUnusedSamplers(configureMCMC(Rmodel))
    }),
    nimbleAFSS = quote({
        assignAFSS(removeUnusedSamplers(configureMCMC(Rmodel)))
    }),
    nimbleRWB = quote({
        assignRWB(removeUnusedSamplers(configureMCMC(Rmodel)))
    })
)

## NIMBLE: We choose to monitor the actual parameters, not derived quantities, for later comparison purposes.
params <- c("theta", "ltheta", "phi", "beta0", "beta", "alpha0", "mean.p", "alpha")
win.data1.compare <- win.data1
win.data1.compare[['nsite']] <- NULL
win.data1.compare[['nrep']] <- NULL
# Compare MCMC efficiency of NIMBLE to jags:

Section6p11p1_compare <- compareMCMCs(
  modelInfo = list(
    code = Section6p11_code,
    data = win.data1.compare,
    constants = list(nsite = 263, nrep = 3),
    inits = initsValues
  ),
  monitors = params,
  MCMCdefs = MCMCdefs,
  MCMCs = c('jags','nimble', 'nimble2', 'nimbleAFSS', 'nimbleRWB'),
  summary = TRUE,
  burnin = nb,
  niter = ni
)

make_MCMC_comparison_pages(Section6p11p1_compare, modelNames = "Section6p11p1", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section6p11p1.html"))

