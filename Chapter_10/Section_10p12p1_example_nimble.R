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
# 10.12.1 Occupancy models with "survival model" observation process: 
#         exponential time-to-detection (TTD) model with simulated data
# ------------------------------------------------------------------------
library(AHMbook)
simOccttd(M = 250, mean.psi = 0.4, mean.lambda = 0.3, beta1 = 1, alpha1 = -1, Tmax = 10)

set.seed(1)
data <- simOccttd()
str(data)

# Plot response (not shown)
hist(data$ttd, breaks = 50, col = "grey", main = "Observed distribution of time to detection", xlim = c(0, data$Tmax), xlab = "Measured time to detection")
abline(v = data$Tmax, col = "grey", lwd = 3)

# Bundle data
str( win.data <- list(ttd = data$ttd, d = data$d, covA = data$covA, 
                      covB = data$covB, nobs = data$M, Tmax = data$Tmax) )

library(nimble)
# Define exponential observation model
# "model1.txt" in book code
Section10p12p1_code <- nimbleCode( {

  # Priors
  int.psi ~ dunif(0, 1)               # Intercept occupancy on prob. scale
  beta1 ~ dnorm(0, 0.001)             # Slope coefficient in logit(occupancy)
  int.lambda ~ dgamma(0.0001, 0.0001) # Poisson rate parameter
  alpha1 ~ dnorm(0, 0.001)            # Slope coefficient in log(rate)
  
  # Likelihood
  for (i in 1:nobs){
  # Model for occurrence
     z[i] ~ dbern(psi[i])
     logit(psi[i]) <- logit(int.psi) + beta1 * covB[i]
     
     # Observation model
     # Exponential model for time to detection ignoring censoring
     ttd[i] ~ dexp(lambda[i])
     log(lambda[i]) <- log(int.lambda) + alpha1 * covA[i]
     # Model for censoring due to species absence and ttd>=Tmax
     d[i] ~ dbern(theta[i])
     theta[i] <- z[i] * step(ttd[i] - Tmax) + (1 - z[i])
  }
  # Derived quantities
  n.occ <- sum(z[1:nobs])                  # Number of occupied sites among M
  }
)

# Inits function for some params
# Initialize with z = 1 throughout and 
# all missings due to censoring, rather than non-occurrence
zst <- rep(1, length(win.data$ttd))
ttdst <-rep(win.data$Tmax+1, data$M)
ttdst[win.data$d == 0] <- NA
inits <- function(){list(z =zst, 
                         ttd = ttdst, 
                         int.psi = runif(1), 
                         int.lambda = runif(1))}

# Parameters to estimate
params <- c("int.psi", "beta1", "int.lambda", "alpha1", "n.occ")

# MCMC settings
ni <- 12000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3

# Call nimble from R  and summarize posteriors
out1 <- nimbleMCMC(
  code = Section10p12p1_code,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

## inspect posterior information/output
library(mcmcplots)
colnames(out1)

Section10p12p1_compare <- compareMCMCs(
  modelInfo = list(
    code = Section10p12p1_code,
    data = win.data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = 2*nb,
  niter = 2*ni
)

make_MCMC_comparison_pages(Section10p12p1_compare, modelNames = "Section10p12p1", dir = outputDirectory)

browseURL(Sys.glob(file.path(outputDirectory, "Section10p12p1.html")))

