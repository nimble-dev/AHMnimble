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
## STOPPED HERE
source("Section_11p6_setup.R")
# 11.7.1 The simplest DR community model with data augmentation
# ------------------------------------------------------------------------
# Augment data set (DA part)
nz <- 150                # Number of potential species in superpopulation
M <- nspec + nz          # Size of augmented data set ('superpopulation')
yaug <- cbind(ysum, array(0, dim=c(nsite, nz))) # Add all zero histories

# Bundle and summarize data set
str( win.data <- list(yaug = yaug, nsite = nrow(ysum), nrep = data$sites$nsurvey, M = M, nspec = nspec, nz = nz) )
library(nimble)
# Specify model in BUGS language
# "model9.txt" in book code
Section11p7p1_code <- nimbleCode({
    
    # Priors to describe heterogeneity among species in community
    for(k in 1:M){                  # Loop over all species in augmented list
      lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
      lp[k] ~ dnorm(mu.lp, tau.lp)
    }
    
    # Hyperpriors to describe full community
    omega ~ dunif(0,1)              # Data augmentation or 'occupancy' parameter
    mu.lpsi ~ dnorm(0,0.001)        # Community mean of occupancy (logit)
    mu.lp ~ dnorm(0,0.001)          # Community mean of detection (logit)
    tau.lpsi <- pow(sd.lpsi, -2)
    sd.lpsi ~ dunif(0,5)            # Species heterogeneity in logit(psi)
    tau.lp <- pow(sd.lp, -2)
    sd.lp ~ dunif(0,5)              # Species heterogeneity in logit(p)
    
    # Superpopulation process:this is the 'paramater expansion' part of PX-DA
    for(k in 1:M){
      w[k] ~ dbern(omega)           # Metacommunity membership indicator
    }                               # (or data augmentation variable)
    
    # Ecological model for latent occurrence z (process model)
    for(k in 1:M){
      mu.psi[k] <- w[k] * psi[k]    # species not part of community zeroed out for z
      logit(psi[k]) <- lpsi[k]
      for (i in 1:nsite) {
        z[i,k] ~ dbern(mu.psi[k])
      }
    }
    
    # Observation model for observed detection frequencies
    for(k in 1:M){
      logit(p[k]) <- lp[k]
      for (i in 1:nsite) {
        mu.p[i,k] <- z[i,k] * p[k]  # non-occurring species are zeroed out for p
        yaug[i,k] ~ dbin(mu.p[i,k], nrep[i])
      }
    }
    
    # Derived quantities
    for(k in 1:M){
      Nocc.fs[k] <- sum(z[1:nsite, k])     # Number of occupied sites among the 267
    }
    for (i in 1:nsite) {
      Nsite[i] <- sum(z[i, 1:M])       # Number of occurring species at each site
    }
    n0 <- sum(w[(nspec+1):(nspec+nz)]) # Number of unseen species in metacommunity
    Ntotal <- sum(w[1:M])              # Total metacommunity size (= nspec + n0)
    }
)


# Initial values
wst <- rep(1, nspec+nz)                   # Simply set everybody at 'occurring'
zst <- array(1, dim = c(nsite, nspec+nz)) # ditto for z
inits <- function() list(z = zst, 
                         w = wst, 
                         lpsi = rnorm(n = nspec+nz), 
                         lp = rnorm(n = nspec+nz))

# Parameters monitored
params <- c("mu.lpsi", "sd.lpsi", "mu.lp", "sd.lp", "Ntotal", "omega", "n0")

# MCMC settings
ni <- 22000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3

# Call nimble from R
out9 <- nimbleMCMC(
  code = Section11p7p1_code,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

## inspect posterior distributions/traceplots
library(mcmcplots)
colnames(out9)
params_of_interest <- params
mcmcplot(out9[, params_of_interest])

Section11p7p1_compare <- compareMCMCs(
  modelInfo = list(
    code = Section11p7p1_code,
    data = win.data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = nb,
  niter = 2*ni
)

make_MCMC_comparison_pages(Section11p7p1_compare, modelNames = "Section11p7p1", dir = outputDirectory)

browseURL(Sys.glob(file.path(outputDirectory, "Section11p7p1.html")))

