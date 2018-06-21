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
# 11.7.2 Dorazio-Royle community model with covariates
# ------------------------------------------------------------------------
# Augment data set: choose one of two different priors on Ntotal
nz <- 250                 # Use for vague prior on Ntotal: M = 395
nz <- 215 - nspec         # Use for informative prior on Ntotal: M = 215
Yaug <- array(0, dim=c(nsite, nrep, nspec+nz)) # array with only zeroes
Yaug[,,1:nspec] <- Y      # copy into it the observed data

# Create same NA pattern in augmented species as in the observed species
missings <- is.na(Yaug[,,1]) # e.g., third survey in high-elevation quads
for(k in (nspec+1):(nspec+nz)){
  Yaug[,,k][missings] <- NA
}
attach(data)
# Bundle and summarize data
str(win.data <- list(Y = Yaug, nsite = dim(Y)[1], nrep = dim(Y)[2], nspec = dim(Y)[3], nz = nz, M = nspec + nz, ele = sites$elev, forest = sites$forest, DAT = date, DUR = dur) )


# Specify model in BUGS language
Section11p7p2_code <- nimbleCode({
    
    # Priors
    omega ~ dunif(0,1)
    # Priors for species-specific effects in occupancy and detection
    for(k in 1:M){
      lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)    # Hyperparams describe community
      betalpsi1[k] ~ dnorm(mu.betalpsi1, tau.betalpsi1)
      betalpsi2[k] ~ dnorm(mu.betalpsi2, tau.betalpsi2)
      betalpsi3[k] ~ dnorm(mu.betalpsi3, tau.betalpsi3)
      lp[k] ~ dnorm(mu.lp, tau.lp)
      betalp1[k] ~ dnorm(mu.betalp1, tau.betalp1)
      betalp2[k] ~ dnorm(mu.betalp2, tau.betalp2)
      betalp3[k] ~ dnorm(mu.betalp3, tau.betalp3)
    }
    
    # Hyperpriors
    # For the model of occupancy
    mu.lpsi ~ dnorm(0,0.01)
    tau.lpsi <- pow(sd.lpsi, -2)
    sd.lpsi ~ dunif(0,8)   # as always, bounds of uniform chosen by trial and error
    mu.betalpsi1 ~ dnorm(0,0.1)
    tau.betalpsi1 <- pow(sd.betalpsi1, -2)
    sd.betalpsi1 ~ dunif(0, 4)
    mu.betalpsi2 ~ dnorm(0,0.1)
    tau.betalpsi2 <- pow(sd.betalpsi2, -2)
    sd.betalpsi2 ~ dunif(0,2)
    mu.betalpsi3 ~ dnorm(0,0.1)
    tau.betalpsi3 <- pow(sd.betalpsi3, -2)
    sd.betalpsi3 ~ dunif(0,2)
    
    # For the model of detection
    mu.lp ~ dnorm(0,0.1)
    tau.lp <- pow(sd.lp, -2)
    sd.lp ~ dunif(0, 2)
    mu.betalp1 ~ dnorm(0,0.1)
    tau.betalp1 <- pow(sd.betalp1, -2)
    sd.betalp1 ~ dunif(0,1)
    mu.betalp2 ~ dnorm(0,0.1)
    tau.betalp2 <- pow(sd.betalp2, -2)
    sd.betalp2 ~ dunif(0,1)
    mu.betalp3 ~ dnorm(0,0.1)
    tau.betalp3 <- pow(sd.betalp3, -2)
    sd.betalp3 ~ dunif(0,1)
    
    # Superpopulation process: Ntotal species sampled out of M available
    for(k in 1:M){
      w[k] ~ dbern(omega)
    }
    
    # Ecological model for true occurrence (process model)
    for(k in 1:M){
      for (i in 1:nsite) {
        logit(psi[i,k]) <- lpsi[k] + betalpsi1[k] * ele[i] + 
        betalpsi2[k] * pow(ele[i],2) + betalpsi3[k] * forest[i]
        mu.psi[i,k] <- w[k] * psi[i,k]
        z[i,k] ~ dbern(mu.psi[i,k])
      }
    }
    
    # Observation model for replicated detection/nondetection observations
    for(k in 1:M){
      for (i in 1:nsite){
        for(j in 1:nrep){
          logit(p[i,j,k]) <- lp[k] + betalp1[k] * DAT[i,j] + 
          betalp2[k] * pow(DAT[i,j],2) + betalp3[k] * DUR[i,j]
          mu.p[i,j,k] <- z[i,k] * p[i,j,k]
          Y[i,j,k] ~ dbern(mu.p[i,j,k])
        }
      }
    }
    
    # Derived quantities
    #for(k in 1:M){
    #   Nocc.fs[k] <- sum(z[,k])       # Number of occupied sites among the 267
    #}
    for (i in 1:nsite){
      Nsite[i] <- sum(z[i,1:M])          # Number of occurring species at each site
    }
    n0 <- sum(w[(nspec+1):(nspec+nz)]) # Number of unseen species
    Ntotal <- sum(w[1:M])                 # Total metacommunity size
    
    # Vectors to save (S for �save�; discard posterior samples for 
    # all minus 1 of the potential species to save disk space)
    # we do this for nz = 250 (i.e., M = 395)
    lpsiS[1:(nspec+1)] <- lpsi[1:(nspec+1)]
    betalpsi1S[1:(nspec+1)] <- betalpsi1[1:(nspec+1)]
    betalpsi2S[1:(nspec+1)] <- betalpsi2[1:(nspec+1)]
    betalpsi3S[1:(nspec+1)] <- betalpsi3[1:(nspec+1)]
    lpS[1:(nspec+1)] <- lp[1:(nspec+1)]
    betalp1S[1:(nspec+1)] <- betalp1[1:(nspec+1)]
    betalp2S[1:(nspec+1)] <- betalp2[1:(nspec+1)]
    betalp3S[1:(nspec+1)] <- betalp3[1:(nspec+1)]
})
    

# Initial values
wst <- rep(1, nspec+nz)                   # Simply set everybody at occurring
zst <- array(1, dim = c(nsite, nspec+nz)) # ditto
inits <- function() list(z = zst, 
                         w = wst, 
                         lpsi = rnorm(n = nspec+nz), 
                         betalpsi1 = rnorm(n = nspec+nz), 
                         betalpsi2 = rnorm(n = nspec+nz), 
                         betalpsi3 = rnorm(n = nspec+nz), 
                         lp = rnorm(n = nspec+nz), 
                         betalp1 = rnorm(n = nspec+nz), 
                         betalp2 = rnorm(n = nspec+nz), 
                         betalp3 = rnorm(n = nspec+nz))
initsValues <- inits()

# Set 1
params1 <- c("omega", "mu.lpsi", "sd.lpsi", "mu.betalpsi1", "sd.betalpsi1", "mu.betalpsi2", "sd.betalpsi2", "mu.betalpsi3", "sd.betalpsi3", "mu.lp", "sd.lp", "mu.betalp1", "sd.betalp1", "mu.betalp2", "sd.betalp2", "mu.betalp3", "sd.betalp3", "Ntotal", "Nsite")

# MCMC settings
ni <- 15000   ;   nt <- 10   ;   nb <- 5000   ;   nc <- 3

# Run nimble, check convergence and summarize posteriors
out101 <- nimbleMCMC(
  code = Section11p7p2_code,
  constants = win.data,
  inits = initsValues,
  monitors = params1,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)
  
