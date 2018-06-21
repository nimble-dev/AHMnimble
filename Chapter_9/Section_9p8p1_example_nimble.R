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
source("Chapter_9_setup.R")
# 9.8.1 Distance sampling with location of encounter
# ------------------------------------------------------------------------

# Simulate a data set and harvest the output
library(AHMbook)
set.seed(1234)
str(tmp <- sim.pdata(N=200,sigma=1,keep.all=FALSE,B=3))

# Harvest some data objects
B <- tmp$B
d <- tmp$d
u1 <- tmp$u1
u2 <- tmp$u2
nind <- length(d)

# Data augmentation
M <- 400                       # max of 400 individuals
nz <- M-nind                   # augment by nz individuals
y <- c(rep(1,nind), rep(0,nz)) # augmented data augmentation variable
u <- cbind(u1,u2)               # Locations of individuals detected
u <- rbind(u, matrix(NA, nrow=nz, ncol=2))


# Bundle and summarize the data set
str(data <- list (B=B, nind=nind, u=u, y=y, nz=nz))

library(nimble)
# Write out the BUGS model file
# "model1.txt" in book code
Section9p8p1_code <- nimbleCode({ 

  # Priors
  sigma ~ dunif(0,10)
  psi ~ dunif(0,1)
   
  # Categorical observation model
  for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)
    u[i,1] ~ dunif(0, 2*B)  # Here is the uniformity assumption made explicit
    u[i,2] ~ dunif(0, 2*B)
    # Compute distance as a derived quantity
    d[i] <- pow( pow( u[i,1]-B,2) + pow(u[i,2]-B,2), 0.5) # Pythagoras
    p[i] <- exp(-d[i]*d[i] / (2*sigma*sigma))
    mu[i] <- p[i] * z[i]
    y[i] ~ dbern(mu[i])
  }
  # Other derived quantity
  N <- sum(z[1:(nind+nz)])
  D <- N / (B*B)
  }
)


# Specify MCMC settings
ni <- 22000   ;   nb <- 2000   ;   nthin <- 2   ;   nc <- 3

# Inits and parameters
inits <- function(){
  list(sigma=runif(1,1,10), 
       psi=runif(1),
       z = c(rep(1,nind), rep(0,nz)) ) }
params <- c("sigma", "N", "psi")

# Execute nimble and summarize the posterior distributions
out1 <- nimbleMCMC(
  code = Section9p8p1_code,
  constants = data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

## look at posterior information
library(mcmcplots)
colnames(out1)
mcmcplot(out1)

Section9p8p1_compare <- compareMCMCs(
  modelInfo = list(
    code = Section9p8p1_code,
    data = data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = nb,
  niter = ni
)

make_MCMC_comparison_pages(Section9p8p1_compare, modelNames = "Section9p8p1", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section9p8p1.html"))
