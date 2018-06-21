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
## DOES NOT WORK
source("Chapter_9_setup.R")
# 9.6 Open HDS models: Implicit Dynamics
# ------------------------------------------------------------------------
library(AHMbook)
# Obtain a data set
set.seed(1236)
str(tmp <- simHDSopen("point", nreps=7, nyears=5, nsites=100, beta.trend=0.2) )
attach(tmp)

#apply(tmp$M.true,2,sum)  # True population size per year

# Define distance class information
delta <- 0.5
nD <- B%/%delta                 # Number of distance classes
midpt <- seq(delta/2, B, delta) # Mid-point of distance intervals

## THIS IS THE SECTINO OF AHM BOOK CODE THAT DOES NOT WORK
# Create the 4-d array
y4d <- array(0, dim=c(nsites, nD, K, nyears))
for(yr in 1:nyears){
  for(rep in 1:K){
    data <- tmp$data[[yr]][[rep]]
    site <- data[,1]
    dclass <- data[,"d"]%/%delta + 1
    ndclass <- B%/%delta
    dclass <- factor(dclass, levels=  1:ndclass)
    y4d[1:nsites,1:nD,rep,yr] <- table(site, dclass)
  }
}


# Bundle and summarize the data set
nobs <- apply(y4d, c(1,3,4), sum)  # Total detections per site and occasion
str( data <- list(y4d=y4d, nsites=nsites, K=K, nD=nD, midpt=midpt, delta=delta, habitat=habitat, B=B, nobs = nobs, T=tmp$nyears) )

library(nimble)
# Define model in BUGS
# "tempemig4d.txt" in book code
Section7p6_code <- nimbleCode( {
    
    # Prior distributions
    beta0 ~ dnorm(0, 0.01)  # Intercept for log(lambda)
    mean.lam <- exp(beta0)
    beta1 ~ dnorm(0, 0.01)  # Coefficient on habitat
    phi ~ dunif(0,1)        # Probability of availability
    sigma ~ dunif(0,5)      # Detection function parameter
    beta.trend ~ dnorm(0, 0.01)
    
    # Construct the multinomial cell probabilities
    for(b in 1:nD){
      log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma) # half-normal 
      f[b] <- (2*midpt[b]*delta)/(B*B)                # radial density function
      cellprobs[b] <- g[b]*f[b]
      cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
    }
    cellprobs[nD+1] <- 1-sum(cellprobs[1:nD])
    for (s in 1:nsites) {
      for (k in 1:K) {
        pdet[s,k] <- sum(cellprobs[1:nD]) # Distance class probabilities
        pmarg[s,k] <- pdet[s,k]*phi       # Marginal probability
      }
    }
    
    for(t in 1:T){                        # Years
      for (s in 1:nsites) {               # Sites
        for (k in 1:K) {                  # Replicates
        # Model part 4: distance class frequencies
        y4d[s,1:nD,k,t] ~ dmulti(cellprobs.cond[1:nD], nobs[s,k,t])
        # Model part 3: total number of detections:
        nobs[s,k,t] ~ dbin(pmarg[s,k], M[s,t])  
        # Model part 2: Availability. Not used in this model but simulated.
        Navail[s,k,t] ~ dbin(phi, M[s,t]) 
       }  # end k loop
      # Model part 1: Abundance model
      M[s,t] ~ dpois(lambda[s,t])
      log(lambda[s,t]) <- beta0 + beta1*habitat[s] + beta.trend*(t-2.5)
      }  # end s loop
    } # end t loop
    
    # Derived quantities
    for(t in 1:T){
      Mtot[t] <- sum(M[1:nsites,t])
    for(k in 1:K){ 
      Ntot[k,t] <- sum(Navail[1:nsites,k,t])
      }
    }
  } # End model
)

# Inits and parameters to save
Navail.st <- apply(y4d, c(1,3,4),sum)
Mst <- apply(Navail.st, c( 1,3), max) +2
inits <- function(){
  list(M=Mst, Navail = Navail.st, sigma = 1.0, phi=.9,beta0=log(2),beta1=.5)
}

params <- c("sigma", "phi", "beta0", "mean.lam", "beta.trend",
   "beta1", "Mtot", "Ntot")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 5   ;   nc <- 3

# Run NIMBLE, look at trace plots and posteriors
outRD <- nimbleMCMC(
  code = Section7p6_code,
  constants = data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

library(coda)
library(mcmcplots)
colnames(outRD)
mcmcplots(outRD)

Section7p6_compare <- compareMCMC(
  modelInfo = list(
    code = Section7p6_code,
    data = data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = nb,
  niter = ni
)

make_MCMC_comparison_pages(Section7p6_compare, modelNames = "Section7p6", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section7p6.html"))

