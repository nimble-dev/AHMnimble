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
# 10.14 Modeling wiggly covariate relationships: penalized splines in hierarchical models
# ------------------------------------------------------------------------
library(AHMbook)
# Execute the function and inspect file produced
data <- wigglyOcc(seed = 1)
str(data)

# Convert matrix data into vectors and prepare stuff
y <- c(data$y)                      # Detection/nondetection data (response)
Xsite <- data$Xsite                 # Fine as is
Xsurvey <- c(data$Xsurvey)          # Survey covariate
site <- rep(1:data$M, data$J)       # Site index

# tmp1 <-spline.prep(Xsite, 20) # This would choose 20 knots for Xsite
tmp1 <-spline.prep(Xsite, NA)   # Choose variable default number of knots
tmp2 <-spline.prep(Xsurvey, NA)
Xocc <- tmp1$X         # Fixed-effects part of covariate Xsite in occ
Zocc <- tmp1$Z         # Random-effects part of covariate Xsite in occ
Xdet <- tmp2$X         # Fixed-effects part of covariate Xsite in det
Zdet <- tmp2$Z         # Random-effects part of covariate Xsite in det
nk.occ <- length(tmp1$knots)    # Number of knots in occupancy spline
nk.det <- length(tmp2$knots)    # Number of knots in detection spline


# Bundle and summarize data set
win.data <- list(y1 = y, 
                 site = site, 
                 M = data$M, 
                 Xocc = Xocc, 
                 Zocc = Zocc, 
                 nk.occ = nk.occ, 
                 Xdet = Xdet, 
                 Zdet = Zdet, 
                 nk.det = nk.det, 
                 nobs = data$M*data$J, 
                 y2 = y, 
                 onesSite = rep(1, 240), 
                 onesSurvey = rep(1, 720), 
                 Xsite = Xsite, 
                 Xsite2 = Xsite^2,  
                 Xsurvey = Xsurvey, 
                 Xsurvey2 = Xsurvey^2)  

str(win.data)    # onesSite and onesSurvey are for occ and det intercepts

library(nimble)
# Specify two models in one in BUGS language
# "hypermodel.txt" in book code
Section10p14_code <- nimbleCode( {
    
    # *** Spline model for the data***
    # --------------------------------
    # Priors
    for(k in 1:3){                 # Regression coefficients
      alpha1[k] ~ dnorm(0, 0.1)   # Detection model
      beta1[k] ~ dnorm(0, 0.1)    # Occupancy model
    }
    for(k in 1:nk.occ){ # Random effects at specified knots (occupancy)
      b.occ[k] ~ dnorm(0, tau.b.occ)
    }
    for(k in 1:nk.det){ # Random effects at specified knots (detection)
      b.det[k] ~ dnorm(0, tau.b.det)
    }
    tau.b.occ ~ dgamma(0.01, 0.01)
    tau.b.det ~ dgamma(0.01, 0.01)
    
    # Likelihood
    # Model for latent occupancy state
    for (i in 1:M) {
      z1[i] ~ dbern(psi1[i])
      logit(psi1[i]) <- fix.terms.occ[i] + smooth.terms.occ[i]
      fix.terms.occ[i] <- beta1[1]*Xocc[i,1] + beta1[2]*Xocc[i,2] + beta1[3]*Xocc[i,2]
      smooth.terms.occ[i] <- inprod(b.occ[1:nk.occ], Zocc[i,])
    }
    # Model for observations
    for(i in 1:nobs){
      y1[i] ~ dbern(mu.y1[i])
      mu.y1[i] <- z1[site[i]] * p1[i]
      logit(p1[i]) <- fix.terms.det[i] + smooth.terms.det[i]
      fix.terms.det[i] <- alpha1[1]*Xdet[i,1] + alpha1[2]*Xdet[i,2] + 
      alpha1[3]*Xdet[i,2]
      smooth.terms.det[i] <- inprod(b.det[1:nk.det], Zdet[i,])
    }
    
    # Derived quantities
    sum.z1 <- sum(z1[1:M])            # Number of occupied sites in sample
    sd.b.occ <- sqrt(1/tau.b.occ)  # SD of spline random effects variance Occ.
    sd.b.det <- sqrt(1/tau.b.det)  # SD of spline random effects variance Det.
    
    
    # *** Polynomial model for same data ***
    # --------------------------------------
    # Priors
    for(k in 1:3){                 # Regression coefficients
      alpha2[k] ~ dnorm(0, 0.1)  # Detection model
      beta2[k] ~ dnorm(0, 0.1)   # Occupancy model
    }
    
    # Likelihood
    # Model for latent occupancy state
    for (i in 1:M) {
    z2[i] ~ dbern(psi2[i])
    logit(psi2[i]) <- beta2[1]*onesSite[i] + beta2[2]*Xsite[i] + beta2[3]*Xsite2[i]
    }
    # Model for observations
    for(i in 1:nobs){
      y2[i] ~ dbern(mu.y2[i])
      mu.y2[i] <- z2[site[i]] * p2[i]
      logit(p2[i]) <- alpha2[1]*onesSurvey[i] + alpha2[2]*Xsurvey[i] + 
      alpha2[3] * Xsurvey2[i]
    }
    
    # Derived quantities
    sum.z2 <- sum(z2[1:M])            # Number of occupied sites in sample
    }
)

# Initial values
zst <- apply(data$y, 1, max)
inits <- function(){list(z1=zst, 
                         alpha1=rnorm(3), 
                         beta1=rnorm(3), 
                         b.occ = rnorm(nk.occ), 
                         b.det = rnorm(nk.det), 
                         tau.b.occ = runif(1), 
                         tau.b.det = runif(1), 
                         z2=zst, 
                         alpha2=rnorm(3), 
                         beta2=rnorm(3))}

# Parameters monitored
params <- c("alpha1", "beta1", "psi1", "p1", "fix.terms.occ", "smooth.terms.occ", "b.occ", "fix.terms.det", "smooth.terms.det", "b.det", "sum.z1", "sd.b.occ", "sd.b.det", "alpha2", "beta2", "psi2", "p2", "sum.z2")

# MCMC settings
ni <- 100000   ;   nb <- 10000   ;   nt <- 90   ;   nc <- 3

# Call nimble from R and summarize posteriors
fhm <- nimbleMCMC(
  code = Section10p14_code,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  thin = nt,
  samplesAsCodaMCMC = TRUE
)

library(mcmcplots)
colnames(fhm)
params_of_interest <- c('alpha1[1]', 'alpha2[1]', 'b.det[1]', 'b.occ[1]', 'beta1[1]', 'beta2[2]')
mcmcplot(fhm[, params_of_interest])  ## inspect posterior summaries
## mixing looks okay

params_for_comparison <- c('alpha1','beta1','b.det','b.occ','tau.b.occ','tau.b.det','alpha2','beta2')

## compare efficiencies between nimble and jags
Section10p14_compare <- compareMCMCs(
  modelInfo = list(
    code = Section10p14_code,
    data = win.data,
    inits = inits()
  ),
  monitors = params_for_comparison,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = 2*nb,
  niter = 2*ni,
  thin = 2
)

make_MCMC_comparison_pages(Section10p14_compare, modelNames = "Section10p14", dir = outputDirectory)

browseURL(Sys.glob(file.path(outputDirectory, "Section10p14.html")))
