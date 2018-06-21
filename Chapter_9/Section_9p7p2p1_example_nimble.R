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
source("Section_9p7_setup.R")
## First we need to simulate data necessary for section 9.7

# 9.7.1 Simulating the ISSJ data over multiple years
# ------------------------------------------------------------------------
# We load the ISSJ data analyzed in chapter 8, package into an unmarked frame  

# 9.7.2.1 The independence model
# ------------------------------------------------------------------------
library(nimble)
# Write out the BUGS model file
# "Sollmann1.txt" in the book code
Section9p7p2p1_code <- nimbleCode({
    
    # Prior distributions
    # Regression parameters
    alpha0 ~ dunif(0,20)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-20,20)
    beta1 ~ dunif(-20,20)
    beta2 ~ dunif(-20,20)
    beta3 ~ dunif(-20,20)
    beta4 ~ dunif(-20,20) # Population trend parameter
    r ~ dunif(0,5)        # NegBin dispersion parameter
    rout <- log(r)
    
    # 'Likelihood'
    for (s in 1:nsites){
      # Linear model for detection function scale
      log(sigma[s]) <- alpha0+alpha1*chap[s]
      # Compute detection probability
      for(k in 1:nD){
        pi[k,s] <- (2*midpt[k]*delta )/(B*B) 
        log(p[k,s]) <- -midpt[k]*midpt[k]/(2*sigma[s]*sigma[s])
        f[k,s] <- p[k,s]*pi[k,s]
        fc[k,s] <- f[k,s]/pcap[s]
        fct[k,s] <- fc[k,s]/sum(fc[1:nD,s])
      }
      pcap[s]<-sum(f[1:nD,s])  # Overall detection probability
      
      # Process model
      for (t in 1:T){
        log(lambda[s,t]) <- beta0 + beta1*chap[s] + beta2*chap2[s] + beta3*elev[s] + beta4*(t - t/2)  # Note trend parameter here
        y[s,t] ~ dbin(pcap[s], N[s,t])
        N[s,t] ~ dnegbin(prob[s,t], r)
        prob[s,t] <- r/(r+lambda[s,t])
      } # End loop over years
    } # End loop over sites
    
    # Distance sampling observation model for observed (binned) distance data
    for(i in 1:nind){
      dclass[i] ~ dcat(fct[1:nD,site[i]]) 
    }
    # Derived parameters
    for(t in 1:6){
      Ntot[t] <- sum(N[1:nsites,t])
      D[t] <- Ntot[t] / (28.27*nsites)   # 300 m point = 28.27 ha
    }
  }
)

# Set up initial values, parameters vector and MCMC settings
Nst <- y+1  # this is for trend model
inits <- function(){list(N=Nst, 
                         beta0=runif(1), 
                         beta1=runif(1), 
                         beta2=runif(1), 
                         beta3=runif(1), 
                         beta4=runif(1), 
                         alpha0=runif(1,3,5), 
                         alpha1=runif(1), 
                         r = 1)}

params <-c('beta0', 'beta1', 'beta2', 'beta3', 'beta4', 'alpha0', 'alpha1', 'Ntot', 'D', 'r')

ni <- 22000   ;   nb <- 2000   ;   nt <- 1   ;   nc <- 3

# Execute NIMBLE, look at convergence and summarize the results
open1 <- nimbleMCMC(
  code = Section9p7p2p1_code,
  constants = data1,
  inits = inits,
  monitors = params,
  nburnin = nb, 
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

## look at posterior information
library(mcmcplots)
colnames(open1)
mcmcplot(open1)
