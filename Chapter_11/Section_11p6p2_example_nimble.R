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
# This model builds slowly in nimble but then the MCMC runs fine.
source("Section_11p6_setup.R")
# 11.6.2 Community occupancy model with bivariate species-specific random effects
# -----------------------------------------------------
# Bundle and summarize data set
str( win.data <- list(ysum = ysum, 
                      M = nrow(ysum), 
                      J = MHB2014$sites$nsurvey, 
                      nspec = dim(ysum)[2], 
                      R = matrix(c(5,0,0,1), ncol = 2), df = 3) )

library(nimble)
# Specify model in BUGS language
# "model6.txt" in book code
Section11p6p2_code <- nimbleCode( {
    
    # Priors
    for(k in 1:nspec){  # Group lpsi and lp together in array eta
      lpsi[k] <- eta[k,1]
      lp[k] <- eta[k,2]
      eta[k, 1:2] ~ dmnorm(mu.eta[1:2], Omega[1:2,1:2])
    }
    # Hyperpriors
    # Priors for mu.lpsi=mu.eta[1] and mu.lp=mu.eta[2]
    # probs = community means of occupancy and detection probability
    for(v in 1:2){
      mu.eta[v] <- log(probs[v] / (1-probs[v]))
      probs[v] ~ dunif(0,1)
    }
    # Prior for variance-covariance matrix
    Omega[1:2, 1:2] ~ dwish(R[1:2,1:2], df)
    Sigma[1:2, 1:2] <- inverse(Omega[1:2,1:2])
    
    # Ecological model for latent occurrence z (process model)
    for(k in 1:nspec){
      logit(psi[k]) <- lpsi[k]   # Must take outside of i loop
      for (i in 1:M) {
        z[i,k] ~ dbern(psi[k])
      }
    }
    
    # Observation model for observed data Y
    for(k in 1:nspec){
      logit(p[k]) <- lp[k]       # Needs to be outside of i loop
      for (i in 1:M) {
        mu.p[i,k] <- z[i,k] * p[k]
        ysum[i,k] ~ dbin(mu.p[i,k], J[i])
     }
    }
    
    # Derived quantities
    rho <- Sigma[1,2] / sqrt(Sigma[1,1] * Sigma[2,2])  # Correlation coefficient
    for(k in 1:nspec){
      Nocc.fs[k] <- sum(z[1:M,k])   # Number of occupied sites among the 267
    }
    for (i in 1:M) {
      Nsite[i] <- sum(z[i,1:nspec])     # Number of occurring species
    }
  }
)



# Initial values
zst <- apply(y, c(1,3), max) # Observed occurrence as starting values for z
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, 
                         Omega = matrix(c(1,0,0,1), ncol = 2), 
                         eta = matrix(0, nrow = win.data$nspec, ncol = 2))

# Parameters monitored
params <- c("mu.eta", "probs", "psi", "p", "Nsite", "Nocc.fs", "Sigma", "rho")

# MCMC settings
ni <- 20000   ;   nt <- 15   ;   nb <- 5000   ;   nc <- 3

m6 <- nimbleModel(
  code = Section11p6p2_code,
  constants = win.data,
  inits = inits(),
  calculate = FALSE) # Faster to skip calculation here

cm6 <- compileNimble(m6)

mcmc6 <- buildMCMC(m6, monitors = params, thin = nt) # Conjugacy is needed for Wishart

cmcmc6 <- compileNimble(mcmc6, project = m6)

out6 <- runMCMC(cmcmc6, niter = ni, samplesAsCodaMCMC = TRUE)

library(mcmcplots)
colnames(out6)
## We may be more interested in N.site[i] or other variables,
## but they would be a lot to plot.
params_of_interest <- m6$expandNodeNames(c("mu.eta","rho", "Sigma"), returnScalarComponents = TRUE)
mcmcplot(out6[, params_of_interest])

Section11p6p2_compare <- compareMCMCs(
  modelInfo = list(
    code = Section11p6p2_code,
    data = win.data,
    inits = inits()
  ),
  monitors = params_of_interest,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = nb,
  niter = 4*ni
) 

make_MCMC_comparison_pages(Section11p6p2_compare, modelNames = "Section11p6p2", dir = outputDirectory)

browseURL(Sys.glob(file.path(outputDirectory, "Section11p6p2.html")))
