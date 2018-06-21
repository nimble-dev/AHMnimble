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
# 10.3 Simulation and analysis of the simplest possible site-occupancy model
# --------------------------------------------------------------------------

# Choose sample sizes and prepare observed data array y
set.seed(24)                  # So we all get same data set
M <- 100                      # Number of sites
J <- 2                        # Number of presence/absence measurements
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Parameter values
psi <- 0.8                    # Probability of occupancy or presence
p <- 0.5                      # Probability of detection

# Generate presence/absence data (the truth)
z <- rbinom(n = M, size = 1, prob = psi)  # R has no Bernoulli

# Generate detection/nondetection data (i.e. presence/absence measurements)
for(j in 1:J){
  y[,j] <- rbinom(n = M, size = 1, prob = z*p)
}

# Look at data
sum(z)                        # True number of occupied sites
sum(apply(y, 1, max))         # Observed number of occupied sites
head(cbind(z=z, y))           # Truth and measurements for first 6 sites

# Bundle data and summarize data bundle
str( win.data <- list(y = y, M = nrow(y), J = ncol(y)) )

library(nimble)
# Specify model in BUGS language
# "model.txt" in book code.
Section10p3_code <- nimbleCode({
   # Priors
   psi ~ dunif(0, 1)
   p ~ dunif(0, 1)
   # Likelihood
   for (i in 1:M) {    # Loop over sites
      z[i] ~ dbern(psi)         # State model
      for (j in 1:J) { # Loop over replicate surveys
         y[i,j] ~ dbern(z[i]*p)  # Observation model 
      }
   }
})


# Initial values
zst <- apply(y, 1, max)       # Avoid data/model/inits conflict
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("psi", "p")

# MCMC settings
ni <- 5000   ;   nt <- 1   ;   nb <- 1000   ;   nc <- 3

# Call nimble and summarize posteriors
library(nimble)
fm2 <- nimbleMCMC(
  code = Section10p3_code,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

library(mcmcplots)
colnames(fm2)
mcmcplot(fm2)
