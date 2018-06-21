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
# This model builds slowly in nimble but then runs fine.
source("Section_11p6_setup.R")

# Bundle and summarize data set
str( win.data <- list(ysum = ysum, M = nrow(ysum), J = MHB2014$sites$nsurvey[1:nsite], nspec = dim(ysum)[2]) )

library(nimble)
# Specify model in BUGS language
# "model5.txt" in book code
Section11p6p1_code <- nimbleCode( {
    
    # Priors
    for(k in 1:nspec){          # Loop over species
      psi[k] ~ dunif(0, 1)
      p[k] ~ dunif(0, 1)
    }
    
    # Ecological model for latent occurrence z (process model)
    for(k in 1:nspec){          # Loop over species
      for (i in 1:M) {         # Loop over sites
        z[i,k] ~ dbern(psi[k])
      }
    }
    
    # Observation model for observed data Y
    for(k in 1:nspec){          # Loop over species
      for (i in 1:M) {
        mup[i,k] <- z[i,k] * p[k]
        ysum[i,k] ~ dbin(mup[i,k], J[i])
      }
    }
    
    # Derived quantities
    for(k in 1:nspec){          # Loop over species
      Nocc.fs[k] <- sum(z[1:M,k]) # Add up number of occupied sites among the 267
    }
    for (i in 1:M) {            # Loop over sites
      Nsite[i] <- sum(z[i,1:nspec])   # Add up number of occurring species at each site
    }
  }
)


# Initial values
zst <- apply(y, c(1,3), max) # Observed occurrence as inits for z
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, 
                         psi = rep(0.4, nspec), 
                         p = rep(0.4, nspec))
# Parameters monitored
params <- c("psi", "p", "Nsite", "Nocc.fs")

# MCMC settings
ni <- 2500   ;   nt <- 2   ;   nb <- 500   ;   nc <- 3

m5 <- nimbleModel(
  code = Section11p6p1_code,
  constants = win.data,
  inits = inits(),
  calculate = FALSE) # Faster to skip calculation here

cm5 <- compileNimble(m5)

mcmc5 <- buildMCMC(m5, monitors = params, thin = nt, useConjugacy = FALSE) # Faster to skip checks for conjugacy

cmcmc5 <- compileNimble(mcmc5, project = m5)

out5 <- runMCMC(cmcmc5, niter = ni, samplesAsCodaMCMC = TRUE)

# Call nimble from R
out5 <- nimbleMCMC(
  code = Section11p6p1_code,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = 2*ni,
  samplesAsCodaMCMC = TRUE,
  check = FALSE
)


library(mcmcplots)
colnames(out5)
params_of_interest <- c('Nocc.fs[50]', 'Nsite[1]', 'p[1]', 'psi[1]') ## small subsample of parameters
mcmcplot(out5[, params_of_interest])

