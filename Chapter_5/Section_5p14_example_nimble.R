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
# 5.14 Random-effects binomial GLM (binomial GLMM)
# ------------------------------------------------

# Bundle data
win.data <- list(y = y, M = nrow(y), J = ncol(y), elev = elev, forest = forest, elev.forest = elev * forest, wind = wind)
str(win.data)

require(nimble)
# Specify model in BUGS language
# Code for the file "RE.Bernoulli.txt" from original book code
section5p14_code <- nimbleCode( {
    
    # Priors
    mu.alpha0 <- logit(mean.theta)              # Random intercepts
    mean.theta ~ dunif(0,1)
    tau.alpha0 <- pow(sd.alpha0, -2)
    sd.alpha0 ~ dunif(0, 10)
    mu.alpha4 ~ dnorm(0, 0.001)                 # Random slope on wind
    tau.alpha4 <- pow(sd.alpha4, -2)
    sd.alpha4 ~ dunif(0, 10)
    
    for(k in 1:3){
      alpha[k] ~ dnorm(0, 0.001)               # Slopes
    }
    
    # Likelihood
    for (i in 1:M){
      alpha0[i] ~ dnorm(mu.alpha0, tau.alpha0) # Intercept random effects
      re00[i] <- alpha0[i] - mu.alpha0         # same zero-centered
      alpha4[i] ~ dnorm(mu.alpha4, tau.alpha4) # Slope random effects
      re04[i] <- alpha4[i] - mu.alpha4         # same zero-centered
    for(j in 1:J){
      y[i,j] ~ dbern(theta[i,j])
      logit(theta[i,j]) <- alpha0[i] + alpha[1] * elev[i] + alpha[2] * forest[i] + alpha[3] * elev.forest[i] + alpha4[i] * wind[i,j]
    }
  }
})

# inits
# For NIMBLE, initial values for alpha[1:3]
# were added to the function from AHM.
inits <- function() list(alpha0 = rnorm(M),
                         alpha4 = rnorm(M),
                         sd.alpha4 = runif(1),
                         mean.theta = runif(1),
                         mu.alpha4 = rnorm(1),
                         sd.alpha0 = runif(1),
                         alpha = rnorm(3)) # Inits

params <- c("mu.alpha0", "sd.alpha0", "alpha0", "alpha", "mu.alpha4", "sd.alpha4", "alpha4", "re00", "re04")                        # Params
ni <- 30000 ; nt <- 25 ; nb <- 5000 ; nc <- 3                 # MCMC settings

# Call nimble from R:
out9 <- nimbleMCMC(code = section5p14_code,
                   constants = win.data,
                   inits = inits,
                   monitors = params,
                   nburnin = nb,
                   niter = ni)
require(coda)
out9.mcmc <- as.mcmc(out9)
colnames(out9.mcmc)
require(mcmcplots)
params_of_interest <- c("mu.alpha0", "sd.alpha0", "alpha[1]", "alpha[2]", "alpha[3]", "mu.alpha4", "sd.alpha4")
mcmcplot(out9.mcmc[,params_of_interest])
## These are not highly mixed.

Section5p14_compare <- compareMCMCs(
  modelInfo = list(
    code = section5p14_code,
    data = win.data,
    inits = inits()
  ), 
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = nb,
  niter = 2*ni) ## double iterations since mixing was not great.

make_MCMC_comparison_pages(Section5p14_compare, modelNames = "Section5p14", dir = outputDirectory)

browseURL(Sys.glob(file.path(outputDirectory, "Section5p14.html")))
