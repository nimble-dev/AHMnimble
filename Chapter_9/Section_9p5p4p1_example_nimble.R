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
# 9.5.4 Fitting temporary emigration HDS models in BUGS
# ------------------------------------------------------------------------

# 9.5.4.1 Simulating a temporary emigration system
# ------------------------------------------------------------------------
library(AHMbook)
## simHDSopen(type="line", 
##            nsites = 100, 
##            mean.lam = 2, 
##            beta.lam = 0, 
##            mean.sig = 1, 
##            beta.sig = 0, 
##            B = 3, 
##            discard0=TRUE, 
##            nreps=2, 
##            phi=0.7, 
##            nyears=5, 
##            beta.trend = 0)

## Obtain a temporary emigration data set 
set.seed(1234)
str(tmp <- simHDSopen("point", nreps=7, nyears=5, nsites=100)  ) 
attach(tmp)


#apply(tmp$M.true,2,sum)

## Define distance class information

delta <- 0.5
nD <- B%/%delta                 # Number of distance classes
midpt <- seq(delta/2, B, delta) # Mid-point of distance intervals

## THE NEXT BLOCK OF CODE FROM THE AHM BOOK CODE DOES NOT SEEM TO WORK
## Create the 4-d array
y4d <- array(0,dim=c(nsites, nD, K, nyears))
for(yr in 1:nyears){
  for(rep in 1:K){
    data <- tmp$data[[yr]][[rep]]
    site <- data[,1]
    dclass <- data[,"d"]%/%delta + 1
    ndclass <- B%/%delta
    dclass <- factor(dclass, levels= 1:ndclass)
    ok <- try(y4d[1:nsites,1:nD,rep,yr] <- table(site, dclass))
    if(inherits(ok, 'try-error')) stop('got it')
  }
}

y3d <- y4d[,,,1]


## Bundle and summarize the data set
nobs <- apply(y3d, c(1,3), sum)  # Total detections per site and occasion
str( data <- list(y3d=y3d, nsites=nsites, K=K, nD=nD, midpt=midpt, delta=delta, habitat=habitat, B=B, nobs = nobs) )


## Define model in BUGS
Section9p5p4p1_code <- nimbleCode({
    # Prior distributions
    beta0 ~ dnorm(0, 0.01)  # Intercept for log(lambda)
    mean.lam <- exp(beta0)
    beta1 ~ dnorm(0, 0.01)  # Coefficient of lambda on habitat
    phi ~ dunif(0,1)        # Probability of availability
    sigma ~ dunif(0.01,5)   # Distance function parameter
    
    # Detection probs for each distance interval and related things
    for(b in 1:nD){
      log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma) # half-normal 
      f[b] <- (2*midpt[b]*delta)/(B*B)    # radial density function
      cellprobs[b] <- g[b]*f[b]
      cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
    }
    cellprobs[nD+1]<- 1-sum(cellprobs[1:nD])
    
    for (s in 1:nsites) {
      for (k in 1:K) {
        pdet[s,k] <- sum(cellprobs[1:nD])   # Distance class probabilities
        pmarg[s,k] <- pdet[s,k]*phi         # Marginal probability
        
        # Model part 4: distance class frequencies
        y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[1:nD], nobs[s,k])
        # Model part 3: total number of detections:
        nobs[s,k] ~ dbin(pmarg[s,k], M[s])
        # nobs[s,k] ~ dbin(pdet[s,k], Navail[s,k]) # Alternative formulation
        # Model part 2: Availability. Not used in this model but simulated.
        Navail[s,k] ~ dbin(phi, M[s]) 
      }  # end k loop
      # Model part 1: Abundance model
      M[s] ~ dpois(lambda[s])    
      log(lambda[s]) <- beta0 + beta1*habitat[s]
    }  # End s loop
    
    # Derived quantities
    Mtot <- sum(M[1:nsites])
    for(k in 1:K){ 
      Ntot[k]<- sum(Navail[1:nsites,k])
    }
  } # End model
)


# Assemble the initial values and parameters to save for JAGS
Navail.st <- apply(y3d, c(1,3),sum)
Mst <- apply(Navail.st, c( 1), max)  +2
inits <- function(){
  list(M=Mst, sigma = 1.0, phi=0.9, beta0=log(2), beta1=0.5)
}
params <- c("sigma", "phi", "beta0", "mean.lam", "beta1", "Mtot", "Ntot")

# MCMC settings
ni <- 60000   ;   nb <- 10000   ;   nt <- 5   ;   nc <- 3

# Run nimble:
# Some additional init may be needed for nimble.
outTE1 <- nimbleMCMC(
  code = Section9p5p4p1_code,
  constants = data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

## Look at posterior output
library(coda)
library(mcmcplots)
colnames(outTE1)
params_of_interest <- c("sigma", "phi", "beta0", "beta1", "mean.lam", "Mtot")
mcmcplot(outTE1[, params_of_interest])
## mixing looks good

Section9p5p4p1_compare <- compareMCMCs(
  modelInfo = list(
    code = Section9p5p4p1_code,
    data = data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c("nimble", "jags"),
  summary = FALSE,
  burnin = nb,
  niter = ni
)

make_MCMC_comparison_pages(Section9p5p4p1_compare, modelNames = "Section9p5p4p1", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section9p5p4p1.html"))


