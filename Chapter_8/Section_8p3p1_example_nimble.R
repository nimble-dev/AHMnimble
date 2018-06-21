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
source("Chapter_8_setup.R")
# 8.3.1 Bayesian analysis of line transect data
# ------------------------------------------------------------------------
# Get data and do data-augmentation
# Observed distances (meters)

x <- c(71.93, 26.05, 58.47, 92.35, 163.83, 84.52, 163.83, 157.33, 
       22.27, 72.11, 86.99, 50.8, 0, 73.14, 0, 128.56, 163.83, 71.85, 
       30.47, 71.07, 150.96, 68.83, 90, 64.98, 165.69, 38.01, 378.21, 
       78.15, 42.13, 0, 400, 175.39, 30.47, 35.07, 86.04, 31.69, 200, 
       271.89, 26.05, 76.6, 41.04, 200, 86.04, 0, 93.97, 55.13, 10.46, 
       84.52, 0, 77.65, 0, 96.42, 0, 64.28, 187.94, 0, 160.7, 150.45, 
       63.6, 193.19, 106.07, 114.91, 143.39, 128.56, 245.75, 123.13, 
       123.13, 153.21, 143.39, 34.2, 96.42, 259.81, 8.72)

B <- 500 # Strip half-width. Larger than max distance
nind <- length(x)

# Analysis of continuous data using data augmentation (DA)
nz <- 200 # Augment observed data with nz = 200 zeroes
y <- c(rep(1, nind), rep(0, nz)) # Augmented inds. have y=0 by definition
x <- c(x, rep(NA, nz)) # Value of distance are missing for the augmented

# Bundle and summarize data set
str( win.data <- list(nind=nind, nz=nz, x=x, y=y, B=B) )

require(nimble)
# Save text file with BUGS model
# This code is called "model1.txt" in the book code.
Section8p3p1_code_model1 <- nimbleCode( {
    
    # Priors
    sigma ~ dunif(0,1000)  # Half-normal scale
    psi ~ dunif(0,1)       # DA parameter
    
    # Likelihood
    for(i in 1:(nind+nz)){
      # Process model
      z[i] ~ dbern(psi)   # DA variables
      x[i] ~ dunif(0, B)  # Distribution of distances
      # Observation model
      logp[i] <- -((x[i]*x[i])/(2*sigma*sigma)) # Half-normal detection fct.
      p[i] <- exp(logp[i])
      mu[i] <- z[i] * p[i]
      y[i] ~ dbern(mu[i]) # Simple Bernoulli measurement error process
      }
      # Derived quantities
      N <- sum(z[1:(nind + nz)]) # Population size
      D <- N / 60                # Density, with A = 60 km^2 when B = 500
  }
)

# Inits
zst <- y
inits <- function(){ list (psi=runif(1), 
                           z=zst, 
                           sigma=runif(1,40,200)) }

# Params to save
params <- c("N", "sigma", "D")

# Experience the raw power of BUGS and summarize marginal posteriors
out <- nimbleMCMC(code = Section8p3p1_code_model1, 
                  constants = win.data, 
                  inits = inits,
                  monitors = params,
                  samplesAsCodaMCMC = TRUE,
                  nburnin = 1000,
                  niter = 11000)

## look at posterior information
require(coda)
require(mcmcplots)
colnames(out)
mcmcplot(out)
# This is not extremely well mixed.

# compare efficiencies of nimble and jags:
Section8p3p1_compare_model1 <- compareMCMCs(
  modelInfo = list(
    code = Section8p3p1_code_model1,
    data = win.data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = 1000,
  niter = 31000  # Increased to obtain better ESS estimates
)

make_MCMC_comparison_pages(Section8p3p1_compare_model1, modelNames = "Section8p3p1_model1", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section8p3p1_model1.html"))


# Analysis of binned data using data augmentation
delta <- 50                # Width of distance bins
xg <- seq(0, B, delta)     # Make the interval cut points
dclass <- x %/% delta + 1  # Convert distances to distance category
nD <- length(xg) -1        # N intervals = length(xg) if max(x) = B

# Bundle data
# Note data changed to include dclass, nG, bin-width delta and midpt
midpt <- xg[-1] - delta/2  # Interval mid-points
str( win.data <- list (nind=nind, nz=nz, dclass=dclass, y=y, B=B,
                       delta=delta, nD=nD, midpt=midpt) )   # Bundle and summarize

# BUGS model specification
# This code is called "model2.txt" in the book code.
Section8p3p1_code_model2 <- nimbleCode({ 
    # Priors
    psi ~ dunif(0, 1)
    sigma ~ dunif(0, 1000)
    
    # Likelihood
    # construct conditional detection probability and Pr(x) for each bin
    for(g in 1:nD){        # midpt = mid point of each cell
      log(p[g]) <- -midpt[g] * midpt[g] / (2 * sigma * sigma)
      pi[g] <- delta / B  # probability of x in each interval
    }
    
    for(i in 1:(nind+nz)){
      z[i] ~ dbern(psi)             # model for individual covariates
      dclass[i] ~ dcat(pi[1:nD])        # population distribution of distance class
      mu[i] <- z[i] * p[dclass[i]]  # p depends on distance class
      y[i] ~ dbern(mu[i])
    }
    # Derived quantities: Population size and density
    N <- sum(z[1:(nind+nz)])
    D <- N / 60
  }
)

# Inits function
zst <- y # DA variables start at observed value of y
inits <- function(){ list (psi=runif(1), 
                           z=zst, 
                           sigma=runif(1,40,200),
                           dclass = c(
                               rep(as.numeric(NA), 73),
                               sample(1:win.data$nD, win.data$nz, replace = TRUE))
                           ) }

# Parameters to save
params <- c("N", "sigma", "D")

# Unleash WinBUGS and summarize posteriors
out.1 <- nimbleMCMC(code = Section8p3p1_code_model2, 
                    constants = win.data, 
                    inits = inits,
                    monitors = params,
                    nburnin = 1000,
                    niter = 11000,
                    samplesAsCodaMCMC = TRUE)

mcmcplot(out.1)
# This is not extremely well mixed.

# compare efficiencies of nimble and jags:
Section8p3p1_compare_model2 <- compareMCMCs(
  modelInfo = list(
    code = Section8p3p1_code_model2,
    data = win.data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = 1000,
  niter = 31000 # Larger sample size gives more accurate ESS estimates
)

make_MCMC_comparison_pages(Section8p3p1_compare_model2, modelNames = "Section8p3p1_model2", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section8p3p1_model2.html"))

