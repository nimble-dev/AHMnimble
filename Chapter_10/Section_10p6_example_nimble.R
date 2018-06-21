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
# 10.6 A model with lots of covariates: use of R function model.matrix with BUGS
# ------------------------------------------------------------------------------
library(AHMbook)
set.seed(148)
data <- simOcc(time.effects = c(0,0), sd.lp = 0, b = 0, show.plot = FALSE)
str(data)                         # Look at data object

# Create habitat factor
hab <- c(rep("A", 90), rep("B", 90), rep("C", 87)) # must have M = 267 sites

# Bundle and summarize data set
HAB <- as.numeric(as.factor(hab))  # Get numeric habitat factor
str( win.data <- list(y = data$y, M = nrow(data$y), J = ncol(data$y), elev = data$elev, forest = data$forest, wind = data$wind, HAB = HAB) )

library(nimble)
# Specify model in BUGS language
# "modelA.txt" in book code
Section10p6_code <- nimbleCode({
    
    # Priors
    mean.p ~ dunif(0, 1)          # Detection intercept on prob. scale
    alpha0 <- logit(mean.p)       #   same on logit scale
    mean.psi ~ dunif(0, 1)        # Occupancy intercept on prob. scale
    beta0 <- logit(mean.psi)      #   same on logit scale
    
    for(k in 1:2){                # 2 terms in detection model
      alpha[k] ~ dnorm(0, 0.1)   # Covariates on logit(detection)
    }
    
    for(k in 1:11){               # 11 terms in occupancy model
      beta[k] ~ dnorm(0, 0.1)    # Covariates on logit(occupancy)
    }
    
    # Likelihood
    for (i in 1:M) {              # Loop over sites
      z[i] ~ dbern(psi[i])
      logit(psi[i]) <- beta0 +                 # occupancy (psi) intercept
      beta[1] * elev[i] +                    # effect of elev
      beta[2] * forest[i] +                  # effect of forest
      beta[3] * equals(HAB[i],2) +           # effect of habitat 2 (= B)
      beta[4] * equals(HAB[i],3) +           # effect of habitat 3 (= C)
      beta[5] * elev[i] * forest[i] +                     # elev:forest
      beta[6] * elev[i] * equals(HAB[i],2) +              # elev:habB
      beta[7] * elev[i] * equals(HAB[i],3) +              # elev:habC
      beta[8] * forest[i] * equals(HAB[i],2) +            # forest:habB
      beta[9] * forest[i] * equals(HAB[i],3) +            # forest:habC
      beta[10] * elev[i] * forest[i] * equals(HAB[i],2) + # elev:forest:habB
      beta[11] * elev[i] * forest[i] * equals(HAB[i],3)   # elev:forest:habC
    for (j in 1:J) {           # Loop over replicates
      y[i,j] ~ dbern(z[i] * p[i,j])        # WinBUGS would need 'straw man' !
      logit(p[i,j]) <- alpha0 +            # detection (p) intercept
      alpha[1] * elev[i] +              # effect of elevation on p
      alpha[2] * wind[i,j]              # effect of wind on p
    }
  }
})

# Inits
inits <- function(){list(z = apply(data$y, 1, max), 
                         mean.psi = runif(1), 
                         mean.p = runif(1), 
                         alpha = rnorm(2), 
                         beta = rnorm(11))}

# Parameters monitored
params <- c("alpha0", "alpha", "beta0", "beta")

# MCMC settings
ni <- 50000   ;   nt <- 10   ;   nb <- 10000   ;   nc <- 3

# Run nimble look at convergence and summarize posteriors
outA <- nimbleMCMC(
  code = Section10p6_code,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

library(mcmcplots)
colnames(outA)
mcmcplot(outA)


# Create design matrix for occupancy covariates and look at it
occDM <- model.matrix(~ data$elev * data$forest * hab)[,-1] # Drop first col.
head(occDM)             # Look at design matrix
str(occDM)

# Bundle and summarize data set
str( win.data <- list(y = data$y, M = nrow(data$y), J = ncol(data$y), elev = data$elev, wind = data$wind, occDM = occDM) )

# Specify model in BUGS language
# "modelB.txt" in book code
Section10p6_codeB <- nimbleCode({

  # Priors
  mean.p ~ dunif(0, 1)          # Detection intercept on prob. scale
  alpha0 <- logit(mean.p)       #   same on logit scale
  mean.psi ~ dunif(0, 1)        # Occupancy intercept on prob. scale
  beta0 <- logit(mean.psi)      #   same on logit scale
  for(k in 1:2){                # 2 terms in detection model
     alpha[k] ~ dnorm(0, 0.1)   # Covariates on logit(detection)
  }
  for(k in 1:11){               # 11 terms in occupancy model
     beta[k] ~ dnorm(0, 0.1)    # Covariates on logit(occupancy)
  }
  
  # Likelihood
  for (i in 1:M) {
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- beta0 + inprod(beta[1:11], occDM[i,])  # slick !
     for (j in 1:J) {
        y[i,j] ~ dbern(z[i] * p[i,j]) # In WinBUGS need 'straw man'
        logit(p[i,j]) <- alpha0 +     # detection (p) intercept
           alpha[1] * elev[i] +       # effect of elevation on p
           alpha[2] * wind[i,j]       # effect of wind on p
     }
  }
})

# Call nimble from R and summarize posteriors
outB <- nimbleMCMC(
  code = Section10p6_codeB,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

colnames(outB)
mcmcplot(outB)

Section10p6_compareB <- compareMCMCs(
  modelInfo = list(
    code = Section10p6_codeB,
    data = win.data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = nb,
  niter = ni
)

make_MCMC_comparison_pages(Section10p6_compareB, modelNames = "Section10p6B", dir = outputDirectory)

browseURL(Sys.glob(file.path(outputDirectory, "Section10p6B.html")))

