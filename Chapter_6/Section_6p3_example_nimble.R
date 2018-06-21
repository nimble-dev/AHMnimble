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
source("Chapter_6_setup_nimble.R")
# 6.3. Simulation and analysis of the simplest possible N-mixture model
# ------------------------------------------------------------------------

# Choose sample sizes and prepare observed data array C
set.seed(24)                # So we all get same data set
M <- 150                    # Number of sites
J <- 2                      # Number of abu. measurements per site (rep. counts)
C <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Parameter values
lambda <- 2.5               # Expected abundance
p <- 0.4                    # Probability of detection (per individual)

# Generate local abundance data (the truth)
N <- rpois(n = M, lambda = lambda)

# Conduct repeated measurements (generate replicated counts)
for(j in 1:J){
  C[,j] <- rbinom(n = M, size = N, prob = p)
}

# Look at data
# The truth ....
table(N)                    # True abundance distribution
sum(N)                      # True total population size at M sites
sum(N>0)                    # True number of occupied sites
mean(N)                     # True mean abundance (estimate of lambda)

# ... and the observations
table(apply(C, 1, max))     # Observed abundance distribution (max count)
sum(apply(C, 1, max))       # Observed total population size at M sites
sum(apply(C, 1, max)>0)     # Observed number of occupied sites
mean(apply(C, 1, max))      # Observed mean "relative abundance"

head(cbind(N=N, count1=C[,1], count2=C[,2])) # First 6 sites

cor(C)[1,2]

# This section using unmarked is not required for running NIMBLE.
# ML analysis of the model using package unmarked:
## library(unmarked)                  # Load package
## umf <- unmarkedFramePCount(y = C)  # Create um data frame
## summary(umf)                       # Summarize
## (fm1 <- pcount(~1 ~1, data = umf)) # Fit model: get estimates on link scale
## backTransform(fm1, "state")        # Get estimates on natural scale
## backTransform(fm1, "det")

# Bayesian analysis of the model using NIMBLE:

require(nimble)
# Bundle and summarize data set
win.data <- list(C = C, M = nrow(C), J = ncol(C))
str(win.data)                      # Look at data

# Specify model in BUGS language:
# This code corresponds to "model1.txt" in the AHM code
Section6p3_code <- nimbleCode( {
 # Priors
   lambda ~ dgamma(0.001, 0.001)
   p ~ dunif(0, 1)
 # Likelihood
   for (i in 1:M) {
      N[i] ~ dpois(lambda)      # State model
      for (j in 1:J) {
         C[i,j] ~ dbin(p, N[i]) # Observation model
      }
   }
 })

# Specify initial values
Nst <- apply(C, 1, max)       # Avoid data/model/inits conflict
inits <- function(){list(N = Nst)}

# Parameters monitored
params <- c("lambda", "p")

# MCMC settings
ni <- 25000   ;   nt <- 20   ;   nb <- 5000   ;   nc <- 3

# Run nimble from R with output as Coda MCMC object
fm2 <- nimbleMCMC(code = Section6p3_code, 
                  constants = win.data, 
                  inits = inits,
                  monitors = params,
                  niter = ni, 
                  nburnin = nb,
                  samplesAsCodaMCMC = TRUE)


require(coda)
require(mcmcplots)
colnames(fm2)
params_of_interest <- c("lambda", "p")
# inspect posterior distributions
mcmcplot(fm2[,params_of_interest])

# Compare MCMC efficiency between NIMBLE and jags
Section6p3_compare <- compareMCMCs(
  modelInfo = list(
    code = Section6p3_code,
    data = win.data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c("nimble", "jags"),
  summary = FALSE,
  burnin = nb,
  niter = ni
)

make_MCMC_comparison_pages(Section6p3_compare, modelNames = "Section6p3", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section6p3.html"))

## Looking at the actual efficiencies:
Section6p3_compare[[1]]$efficiency

