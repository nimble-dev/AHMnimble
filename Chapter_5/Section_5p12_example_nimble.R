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
# 5.12 Moment-matching in a binomial GLM to accommodate underdispersion
# ---------------------------------------------------------------------

# Bundle data
win.data <- list(y = data$C, M = nrow(data$C), J = ncol(data$C), elev = elev, N = 32)

require(nimble)
# Specify model in BUGS language
# Code for "squeezed_count_GLM.txt" from original book code.
Section5p12_code <- nimbleCode( {
    
    # Priors
    alpha ~ dnorm(0, 1.0E-06)
    beta ~ dnorm(0, 1.0E-06)
    
    # Likelihood
    for (i in 1:M){
    mu[i] <- alpha + beta * elev[i] # linear model for expected response
    logit(p[i]) <- logit(mu[i] / N) # express param as function of first moment
    for(j in 1:J){
    y[i,j] ~ dbin(p[i], N)
    }
  }
})

# Initial values
inits <- function() list(alpha = 1.7, beta = -1.2)          # works always
inits.1 <- function() list(alpha = runif(1), beta = runif(1)) # works sometimes

# Parameters monitored
params <- c("alpha", "beta", "mu")

# MCMC settings
ni <- 3000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call nimble from R
out7 <- nimbleMCMC(code = Section5p12_code,
                   constants = win.data,
                   inits = inits,
                   monitors = params,
                   nburnin = nb,
                   niter = ni)
require(coda)
out7.mcmc <- as.mcmc(out7)
require(mcmcplots)
mcmcplot(out7.mcmc[,1:2])

# Compare MCMC efficiency between NIMBLE and jags:
Section5p12_compare <- compareMCMCs(
  modelInfo = list(
    code = Section5p12_code,
    data = win.data,
    inits = inits()
    ),
  MCMCs = c('nimble', 'jags'),
  summary = TRUE,
  burnin = nb,
  niter = ni
)

make_MCMC_comparison_pages(Section5p12_compare, modelNames = "Section5p12", dir = outputDirectory)

browseURL(Sys.glob(file.path(outputDirectory, "Section5p12.html")))
