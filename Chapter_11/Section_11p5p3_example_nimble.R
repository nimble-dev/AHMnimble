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
# 11.5.3 N-mixture model for the observed community size
# -----------------------------------------------------
# Bundle and summarize data set
str( win.data <- list(CC = CC, M = nrow(CC), J = ncol(CC), ele = ele, forest = forest, DAT = DAT, DUR = DUR) )

library(nimble)
# Specify model in BUGS language
# "model4.txt" in book code
Section11p5p3_code <- nimbleCode({
    
    # Priors
    alpha0 ~ dnorm(0, 0.01)      # Base-line community detection probability
    beta0 ~ dnorm(0, 0.01)       # Base-line community size (number of species)
    for(v in 1:3){
      alpha[v] ~ dnorm(0, 0.01) # Covariate effects on detection
      beta[v] ~ dnorm(0, 0.01)  # Covariate effects on community size
    }
    
    # Likelihood
    # Ecological model for true community size
    for (i in 1:M){              # Loop over sites
      N[i] ~ dpois(lambda[i])   # Community size
      lambda[i] <- exp(beta0 + beta[1] * ele[i] + beta[2] * pow(ele[i],2) + 
      beta[3] * forest[i])
      
      # Observation model for repeated measurements
      for (j in 1:J){          # Loop over occasions
        CC[i,j] ~ dbin(p[i,j], N[i])
        p[i,j] <- 1 / (1 + exp(-lp[i,j]))
        lp[i,j] <- alpha0 + alpha[1] * DAT[i,j] + alpha[2] * pow(DAT[i,j],2) +
        alpha[3] * DUR[i,j]
        # logit(p) = ... causes undefined real result in WinBUGS (but not JAGS)
        }
     }
  }
)


# Define function to generate random initial values
Nst <- apply(CC, 1, max, na.rm = TRUE) + 1
Nst[Nst == -Inf] <- max(Nst, na.rm = T)  # Some nonzero val. for unsurv. sites
inits <- function() list(N = Nst, 
                         alpha0 = rnorm(1), 
                         alpha = rnorm(3), 
                         beta0 = rnorm(1), 
                         beta = rnorm(3))

# Parameters monitored
params <- c("alpha0", "alpha", "beta0", "beta")

# MCMC settings
ni <- 6000   ;   nt <- 4   ;   nb <- 2000   ;   nc <- 3

# Run nimble from R, look at traceplots 
#   and summarize posteriors
out4 <- nimbleMCMC(
 code = Section11p5p3_code,
 constants = win.data,
 inits = inits,
 monitors = params,
 nburnin = nb,
 niter = ni,
 samplesAsCodaMCMC = TRUE
)

library(mcmcplots)
colnames(out4)
mcmcplot(out4)

Section11p5p3_compare <- compareMCMCs(
  modelInfo = list(
    code = Section11p5p3_code,
    data = win.data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = nb,
  niter = 4*ni 
)

make_MCMC_comparison_pages(Section11p5p3_compare, modelNames = "Section11p5p3", dir = outputDirectory)

browseURL(Sys.glob(file.path(outputDirectory, "Section11p5p3.html")))
