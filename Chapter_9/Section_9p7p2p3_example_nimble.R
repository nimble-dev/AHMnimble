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
# 9.7.2.3 The Glorious Integrated HDS/Dail-Madsen Model
# ------------------------------------------------------------------------
library(nimble)
# model "Sollmann3.txt" in book code
Section9p7p2p3_code <- nimbleCode({
    # Prior distributions
    # Regression parameters
    alpha0 ~ dunif(0,20)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-20,20)
    beta1 ~ dunif(-20,20)
    beta2 ~ dunif(-20,20)
    beta3 ~ dunif(-20,20)
    
    # Priors for dynamics parameters: here they are constant across years
    # We could add covariate models for logit(phi) and log(gamma)
    phi ~ dunif(0,1)
    gamma ~ dunif(0,5)
    
    # NegBin dispersion parameter
    r ~ dunif(0,5)
    rout <- log(r)
    
    # 'Likelihood'
    for (s in 1:nsites){
      # Linear model for detection function scale
      log(sigma[s]) <- alpha0+alpha1*chap[s]
      
      # Compute detection probability
      for(k in 1:nD){
        log(p[k,s]) <- -midpt[k]*midpt[k]/(2*sigma[s]*sigma[s])
        f[k,s] <- p[k,s]*pi[k,s]
        fc[k,s] <- f[k,s]/pcap[s]
        fct[k,s] <- fc[k,s]/sum(fc[1:nD,s])
        pi[k,s] <- (2*midpt[k]*delta )/(B*B)
      }
      pcap[s]<-sum(f[1:nD,s])  # Overall detection probability
      
      # Process model
      # Abundance model for year 1
      log(lambda[s,1]) <- beta0 + beta1*chap[s] + beta2*chap2[s] + beta3*elev[s]
      y[s,1] ~ dbin(pcap[s], N[s,1])
      N[s,1] ~ dnegbin(prob[s,1], r)
      prob[s,1] <- r/(r+lambda[s,1])
      
      # Population dynamics model for subsequent years
      for (t in 2:T){                      # Loop over years
        S[s,t] ~ dbinom(phi, N[s, t-1])   # Survivors
        R[s,t] ~ dpois(gamma * N[s, t-1]) # Recruits
        N[s,t] <- S[s,t] + R[s,t]         # N = Survivors + Recruits
        y[s,t]~ dbin(pcap[s],N[s,t])       # Measurement error
      }
    }
    
    # Distance sampling observation model for observed (binned) distance data
    for(i in 1:nind){
      dclass[i] ~ dcat(fct[1:nD,site[i]])
    }
    
    # Derived parameters
    for(t in 1:6){
      Ntot[t] <- sum(N[1:nsites,t])
      D[t] <- Ntot[t] / (28.27*nsites)     # 300 m point = 28.27 ha
      } 
    }
)


# Set up some sensible starting values for S and R
yin <- y+1
yin[,2:6] <- NA
Sin <- Rin <- matrix(NA, nrow=nsites, ncol=nyrs)
y1 <- y + 1

for(s in 1:nsites){
  for (t in 2:6){
    Sin[s,t] <- rbinom(1,y1[s,t-1], phi )
    Rin[s,t] <- ifelse((y1[s,t]-Sin[s,t])>0, y1[s,t]-Sin[s,t], 0)
  }
}

# Set up initial values, parameters vector and MCMC settings
inits <-function(){list(N=yin, beta0=runif(1), beta1=runif(1), beta2=runif(1),
                        beta3=runif(1), alpha0=runif(1,3,5), alpha1=runif(1), phi=0.6, gamma=0.3, 
                        R=Rin, S=Sin) }

params <- c('beta0', 'beta1', 'beta2', 'beta3', 'alpha0', 'alpha1', 'phi', 
            'gamma', 'Ntot', 'D', 'r')

ni <- 152000   ;   nb <- 2000   ;   nt <- 10   ;   nc <- 5

# Run nimble, look at convergence and summarize the results
library(nimble)
open3  <- nimbleMCMC(
  code = Section9p7p2p3_code,
  constants = data1,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

library(mcmcplots)
colnames(open3)
mcmcplot(open3)

