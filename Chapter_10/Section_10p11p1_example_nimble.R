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
source("Chapter_10_setup.R")
# 10.11.1 A magical covariate
# ------------------------------------------------------------------------

library(AHMbook)
set.seed(1)
data <- sim3Occ(nunit = 500,
                nsubunit = 5, 
                nrep = 1, 
                mean.psi = 0.8, 
                beta.Xpsi = 1, 
                sd.logit.psi = 0.4, 
                mean.theta = 0.6, 
                theta.time.range = c(0, 0), 
                beta.Xtheta = 1, 
                sd.logit.theta = 0.6, 
                mean.p = 0.4, 
                p.time.range = c(0,0), 
                beta.Xp = -1, 
                sd.logit.p = 0.8)


# Bundle and summarize data set
y <- data$y
str( win.data <- list(y = y, 
                      nunit = dim(y)[1], 
                      nsubunit = dim(y)[2], 
                      nrep = dim(y)[3], 
                      covA = data$covA, 
                      covB = data$covB, 
                      covC = data$covC) )

library(nimble)
# Define model in BUGS langauge
# "model.txt" in original book code
Section10p11p1_code <- nimbleCode( { 
    
    # Priors 
    int.psi ~ dunif(0,1)   # Occupancy probability
    int.theta ~ dunif(0,1) # Availability probability
    int.p ~ dunif(0,1)     # Detection probability
    beta.lpsi ~ dnorm(0, 0.01)
    beta.ltheta ~ dnorm(0, 0.01)
    beta.lp ~ dnorm(0, 0.01)
    
    # Likelihood
    for (i in 1:nunit){ 
      # Occupancy model for quad i
      z[i] ~ dbern(psi[i])
      logit(psi[i]) <- logit(int.psi) + beta.lpsi * covA[i]
      for (j in 1:nsubunit){
        # Availability in subquad j
        a[i,j] ~ dbern(mu.a[i,j])
        mu.a[i,j] <- z[i] * theta[i,j]
        logit(theta[i,j]) <- logit(int.theta) + beta.ltheta * covB[i,j]
        for (k in 1:nrep){
          # PCR detection error process in replicate k
          y[i,j,k] ~ dbern(mu.y[i,j,k])
          mu.y[i,j,k] <- a[i,j] * p[i,j]
          logit(p[i,j]) <- logit(int.p) + beta.lp * covC[i,j,1]
        }
      }
      tmp[i] <- step(sum(a[i,1:nsubunit])-0.1)
    }
    
    # Derived quantities
    sum.z <- sum(z[1:nunit])      # Total number of occupied quadrats
    sum.a <- sum(tmp[1:nunit])    # Total number of quads with presence in samples
    p.theta <- int.p * int.theta # What a 2-level model estimates as 'p'
    }
)


# Initial values
inits <- function() list(z = array(1, dim = data$nunit), 
                         a = array(1, dim =c(data$nunit, data$nsubunit)) )    # Set all to 1 to avoid conflict

# Parameters monitored
params <- c("int.psi", "int.theta", "int.p", "beta.lpsi", "beta.ltheta", "beta.lp","p.theta", "sum.z", "sum.a") 

# MCMC settings
ni <- 25000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3

# Call nimble and summarize posterior
out <- nimbleMCMC(
  code = Section10p11p1_code,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

library(mcmcplots)
colnames(out)
mcmcplot(out)
## mixing needs to be checked

## compare efficiency to jags:
Section10p11p1_compare <- compareMCMCs(
  modelInfo = list(
    code = Section10p11p1_code,
    data = win.data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = 3*nb,
  niter = 3*ni
)

make_MCMC_comparison_pages(Section10p11p1_compare, modelNames = "Section10p11p1", dir = outputDirectory)

browseURL(Sys.glob(file.path(outputDirectory, "Section10p11p1.html")))

