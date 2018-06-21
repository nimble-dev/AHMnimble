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
source("Section_7p6_setup.R")
#  7.6 Bayesian Analysis in BUGS using the Conditional Multinomial (3-part model)
# ------------------------------------------------------------------------

require(nimble)
# Write BUGS model
# File "model.txt" of original book code
Section7p6_code <- nimbleCode( {
    
    # Prior distributions
    p0 ~ dunif(0,1)
    alpha0 <- logit(p0)
    alpha1 ~ dnorm(0, 0.01)
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    beta2 ~ dnorm(0, 0.01)
    beta3 ~ dnorm(0, 0.01)
    
    for(i in 1:M){ # Loop over sites
      # Conditional multinomial cell probabilities
      pi[i,1] <- p[i]
      pi[i,2] <- p[i]*(1-p[i])
      pi[i,3] <- p[i]*(1-p[i])*(1-p[i])
      pi[i,4] <- p[i]*(1-p[i])*(1-p[i])*(1-p[i])
      pi0[i] <- 1 - (pi[i,1] + pi[i,2] + pi[i,3] + pi[i,4])
      pcap[i] <- 1 - pi0[i]
      for(j in 1:4){
        pic[i,j] <- pi[i,j] / pcap[i]
     }
    
      # logit-linear model for detection: understory cover effect
      logit(p[i]) <- alpha0 + alpha1 * X[i,1]
      
      # Model specification, three parts:
      y[i,1:4] ~ dmulti(pic[i,1:4], n[i]) # component 1 uses the conditional
                                          #    cell probabilities
      n[i] ~ dbin(pcap[i], N[i])          # component 2 is a model for the
                                          #    observed sample size 
      N[i] ~ dpois(lambda[i])             # component 3 is the process model
      
      # log-linear model for abundance: UFC + TRBA + UFC:TRBA
      log(lambda[i])<- beta0 + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,2]*X[i,1]
    }
  }
)

# Parameters monitored

params <- c("p0", "alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3")

# MCMC settings
nc <- 3   ;   ni <- 12000   ;   nb <- 1000   ;   nt <- 1

# Call nimble from R 
out <- nimbleMCMC(code = Section7p6_code, 
                  constants = data, 
                  inits = initsValues, 
                  monitors = params,
                  nburnin = nb, 
                  niter = ni,
                  samplesAsCodaMCMC = TRUE)

## look at posterior output
require(coda)
colnames(out)
require(mcmcplots)
params_of_interest <- c("alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3", "p0")
mcmcplot(out[, params_of_interest])

# Compare MCMC efficiencies with nimble and jags
Section7p6_compare <- compareMCMCs(
  modelInfo = list(
    code = Section7p6_code,
    data = data,
    inits = initsValues,
    monitors = params
  ),
  MCMCs = c('nimble', 'jags', 'nimble_slice'),
  summary = FALSE,
  burnin = nb,
  niter = ni
)

make_MCMC_comparison_pages(Section7p6_compare, modelNames = "Section7p6", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section7p6.html"))


