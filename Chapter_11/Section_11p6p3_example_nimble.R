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
# The models here build slowly in nimble.
source("Section_11p6_setup.R")
# 11.6.3 Modeling species-specific effects in community occupancy models
# ------------------------------------------------------------------------
# Look at distribution of body mass among 145 observed species
mass <- tapply(MHB2014$species$body.mass, MHB2014$species$specid, mean) # Get mean species mass (all 158 spec)
mass <- mass[-toss.out]                           # Only retain 145 observed species
hist(log10(mass), breaks = 40, col = "grey")      # Look at log10
gmass <- as.numeric(log10(mass) %/% 1.3 + 1)      # size groups 1, 2 and 3
gmass[gmass == 4] <- 3                            # Mute swan is group 3, too

# Bundle and summarize data set
str( win.data <- list(ysum = ysum, 
                      g = gmass, 
                      M = nrow(ysum), 
                      J = MHB2014$sites$nsurvey[1:nsite], 
                      nspec = dim(ysum)[2]) )

library(nimble)
# Specify model in BUGS language
# "model7.txt" in book code
Section11p6p3_code <- nimbleCode ({
    
    # Priors
    for(k in 1:nspec){      # loop over species
      lpsi[k] ~ dnorm(mu.lpsi[g[k]], tau.lpsi[g[k]]) # note g-dependence now
      lp[k] ~ dnorm(mu.lp[g[k]], tau.lp[g[k]])
    }
    
    # Hyperpriors
    # We replaced g with gid as the loop index, because g is a variable above.
    for(gid in 1:3){          # loop over groups (gid)
      mu.lpsi[gid] <- logit(mu.psi[gid])      # everythingid is indexed g now
      mu.lp[gid] <- logit(mu.p[gid])
      mu.psi[gid] ~ dunif(0,1)
      mu.p[gid] ~ dunif(0,1)
      tau.lpsi[gid] <- pow(sd.lpsi[gid], -2)
      sd.lpsi[gid] ~ dunif(0,5)
      tau.lp[gid] <- pow(sd.lp[gid], -2)
      sd.lp[gid] ~ dunif(0,5)
    }
    
    # Ecological model for latent occurrence z (process model)
    for(k in 1:nspec){      # no change at all down here in model
      logit(psi[k]) <- lpsi[k]
      for (i in 1:M) {
        z[i,k] ~ dbern(psi[k])
      }
    }
    
    # Observation model for observed data ysum
    for(k in 1:nspec){      # Loop over species
      logit(p[k]) <- lp[k]
      for (i in 1:M) {
        mu.px[i,k] <- z[i,k] * p[k]  # call mu.px to avoid conflict with above
        ysum[i,k] ~ dbin(mu.px[i,k], J[i])
      }
    }
    
    # Derived quantities
    for(k in 1:nspec){          # Loop over species
      Nocc.fs[k] <- sum(z[1:M,k]) # Number of occupied sites among the 267
    }
    for (i in 1:M) {            # Loop over sites
      Nsite[i] <- sum(z[i,1:nspec])   # Number of occurring species at each site
    }
  }
)


# Initial values
zst <- apply(y, c(1,3), max)
zst[is.na(zst)] <- 1
inits <- function() list(z = zst)

# Parameters monitored
params <- c("mu.psi", "mu.lpsi", "sd.lpsi", "mu.p", "mu.lp", "sd.lp")

# MCMC settings
ni <- 6000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3

# Call nimble from R, look at convergence and summarize posteriors
out7 <- nimbleMCMC(
  code = Section11p6p3_code,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

library(mcmcplots)
colnames(out7)
mcmcplot(out7)

# Bundle and summarize data set
logmass <- as.numeric(log10(mass))         # Take log10 of body mass
str( win.data <- list(ysum = ysum,
                      logmass = logmass,
                      M = nrow(ysum),
                      nsite = nrow(ysum), 
                      J = MHB2014$sites$nsurvey[1:nsite],
                      nspec = dim(ysum)[2]) )


# Specify model in BUGS language
# "model8.txt" in book code
Section11p6p3_code_m8 <- nimbleCode( {
    
    # Priors
    for(k in 1:nspec){              # loop over species
      lpsi[k] ~ dnorm(mu.lpsi[k], tau.lpsi[k]) # now all indexed by k, not g
      tau.lpsi[k] <- 1/var.lpsi[k]
      lp[k] ~ dnorm(mu.lp[k], tau.lp[k])
      tau.lp[k] <- 1/var.lp[k]
      mu.lpsi[k] <- delta0.lpsi + delta1.lpsi * logmass[k]
      mu.lp[k] <- delta0.lp + delta1.lp * logmass[k]
      log(var.lpsi[k]) <- phi0.lpsi + phi1.lpsi * logmass[k]
      log(var.lp[k]) <- phi0.lp + phi1.lp * logmass[k]
    }
    # Priors for regression params for means
    delta0.lpsi ~ dnorm(0, 0.01)
    delta1.lpsi ~ dnorm(0, 0.01)
    delta0.lp ~ dnorm(0, 0.01)
    delta1.lp ~ dnorm(0, 0.01)
    # Priors for regression params for variances
    phi0.lpsi ~ dnorm(0, 0.01)
    phi1.lpsi ~ dnorm(0, 0.01)
    phi0.lp ~ dnorm(0, 0.01)
    phi1.lp ~ dnorm(0, 0.01)
    
    # Ecological model for latent occurrence z (process model)
    for(k in 1:nspec){
      logit(psi[k]) <- lpsi[k]
      for (i in 1:M) {
        z[i,k] ~ dbern(psi[k])
      }
    }
    
    # Observation model for observed data ysum
    for(k in 1:nspec){              # Loop over species
      logit(p[k]) <- lp[k]
      for (i in 1:M) {
        mu.p[i,k] <- z[i,k] * p[k]
        ysum[i,k] ~ dbin(mu.p[i,k], J[i])
      }
    }
    
    # Derived quantities
    for(k in 1:nspec){          # Loop over species
      Nocc.fs[k] <- sum(z[1:nsite,k]) # Number of occupied sites among the 267
    }
    for (i in 1:nsite) {        # Loop over sites
      Nsite[i] <- sum(z[i,1:nspec])   # Number of occurring species at each site
    }
  }
)


# Initial values
zst <- apply(y, c(1,3), max)
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, delta0.lpsi = rnorm(1), delta1.lpsi = rnorm(1),
                         delta0.lp = rnorm(1), delta1.lp = rnorm(1), phi0.lpsi = rnorm(1),
                         phi1.lpsi = rnorm(1), phi0.lp = rnorm(1), phi1.lp = rnorm(1))

# Parameters monitored
params <- c("delta0.lpsi", "delta1.lpsi", "delta0.lp", "delta1.lp", "phi0.lpsi",
            "phi1.lpsi", "phi0.lp", "phi1.lp", "psi", "p", "Nocc.fs", "Nsite")

# MCMC settings
ni <- 12000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3

# Call NIMBLE from R, look at convergence and summarize posteriors
out8 <- nimbleMCMC(
  code = Section11p6p3_code_m8,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

colnames(out8)
params_of_interest <- c("delta0.lp", "delta0.lpsi", "delta1.lp", "delta1.lpsi", "phi0.lp", "phi0.lpsi", "phi1.lp", "phi1.lpsi")
mcmcplot(out8[, params_of_interest])

Section11p6p3_compare <- compareMCMCs(
  modelInfo = list(
    code = Section11p6p3_code_m8,
    data = win.data,
    inits = inits()
  ),
  monitors = params_of_interest,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = nb,
  niter = 2*ni
)

make_MCMC_comparison_pages(Section11p6p3_compare, modelNames = "Section11p6p3", dir = outputDirectory)

browseURL(Sys.glob(file.path(outputDirectory, "Section11p6p3.html")))

