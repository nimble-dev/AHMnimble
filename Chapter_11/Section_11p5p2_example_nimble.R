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
source("Section_11p5_setup.R")
# 11.5.2 Poisson random effects model for the observed community size
# ------------------------------------------------------------------------

# Bundle and summarize data set
str( win.data <- list(CC = CC, M = nrow(CC), J = ncol(CC), ele = ele, forest = forest, DAT = DAT, DUR = DUR) )

library(nimble)
# Specify model in BUGS language
# "model3.txt" in book code
Section11p5p2_code <- nimbleCode( {

  # Priors
  mugamma0 ~ dnorm(0, 0.001)     # Hyperparameters
  taugamma0 <- pow(sd.gamma0,-2)
  sd.gamma0 ~ dunif(0, 10)
  for(v in 1:6){                 # Parameters
     gamma[v] ~ dnorm(0, 0.001)
  }
  
  # Likelihood for Poisson GLMM
  for (i in 1:M){                # Loop over sites
    gamma0[i] ~ dnorm(mugamma0, taugamma0)     # site intercepts random now
    for(j in 1:J){               # Loop over repeated measurements
      CC[i,j] ~ dpois(lambda[i,j])
      log(lambda[i,j]) <- gamma0[i] + gamma[1]*ele[i] + gamma[2] * pow(ele[i],2) +
        gamma[3] * forest[i] + gamma[4] * DAT[i,j] + gamma[5] * pow(DAT[i,j],2) +
        gamma[6] * DUR[i,j]
     }
   }
 }
)


# Initial values
inits <- function() list(gamma0 = rnorm(nrow(CC)), gamma = rnorm(6))

# Parameters monitored
params <- c("mugamma0", "sd.gamma0", "gamma")

# MCMC settings
ni <- 6000   ;   nt <- 4   ;   nb <- 2000   ;   nc <- 3

# Call nimble from R 
out3 <- nimbleMCMC(
  code = Section11p5p2_code,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

library(mcmcplots)
colnames(out3)
mcmcplot(out3)

Section11p5p2_compare <- compareMCMCs(
  modelInfo = list(
    code = Section11p5p2_code,
    data = win.data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = 4*nb,
  niter = 4*ni
)

make_MCMC_comparison_pages(Section11p5p2_compare, modelNames = "Section11p5p2", dir = outputDirectory)

browseURL(Sys.glob(file.path(outputDirectory, "Section11p5p2.html")))

