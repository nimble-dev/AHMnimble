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
source("Chapter_7_setup.R")
# 7.8.4 Bayesian analysis using data augmentation (DA): heterogeneity models in BUGS
# ----------------------------------------------------------------------------------
# Extract data and do data augmentation up to M = 400

alfl <- read.csv(system.file("csv", "alfl.csv", package="unmarked"))
head(alfl, 5)

alfl.covs <- read.csv(system.file("csv", "alflCovs.csv", package="unmarked"),
                      row.names=1)
head(alfl.covs)

y <- as.matrix(alfl[,c("interval1","interval2","interval3")] )
nind <- nrow(y)
M <- 400
y <- rbind(y, matrix(0, nrow=(M-nind), ncol=3))

# Site ID
# Make site ID into an integer: This only works ok here because sites are in 
#  alphabetical order in the data set!
site <- as.numeric(alfl$id)
site <- c(site, rep(NA, M-nind))

# Next we extract the covariates and standardize them
sitecovs <- scale(alfl.covs[,c("woody", "struct", "time.1")])


# This model contains dynamic indexing of X.
# To make this work in NIMBLE, it is necessary
# to distinguish between constants and data
data <- list(y = y, X = sitecovs, group = site)
str(data)
constants <- list( J = 3, M = 400, nsites = 50 )
str(constants)

require(nimble)
# Specify model in BUGS language
Section7p8p4_code <- nimbleCode({
    
  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0 / (1-p0))    # same as logit(p0)	
  alpha1 ~ dnorm(0, 0.01)
  alpha2 ~ dnorm(0,0.01)
  beta0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)
  psi <- sum(lambda[1:nsites]) / M   # psi is a derived parameter
  
  # log-linear model for abundance: lambda depends on WOODY
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * X[s,1]
    probs[s] <- lambda[s] / sum(lambda[1:nsites])
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[1:nsites])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Observation model: p depends on 2 covariates: STRUCT + TIME  
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * X[group[i],2] + alpha2*X[group[i],3]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
    }
  }
})

# Parameters monitored
params <- c("p0", "alpha0", "alpha1", "alpha2", "beta0", "beta1", "psi")

# Initial values
inits <- function() {
  list (p0 = runif(1), 
        alpha1 = runif(1), 
        alpha2 = rnorm(1),
        beta0=runif(1), 
        beta1=rnorm(1), 
        z= c( rep(1,100), 
             rep(0, 300)),
        # To avoid nimble warnings, we initialize unobserved groups
        group = c(rep(NA,length(alfl$id)),
                  sample(1:nrow(data$X), size = length(site) - length(alfl$id), replace = TRUE))
  )}

# MCMC settings
ni <- 11000   ;   nb <- 1000   ;   nt <- 4   ;   nc <- 3

out <- nimbleMCMC(code = Section7p8p4_code, 
                  data = data, 
                  constants = constants,
                  inits = inits,
                  monitors = params,
                  nburnin = nb, 
                  niter = ni,
                  samplesAsCodaMCMC = TRUE)
require(coda)
require(mcmcplots)
colnames(out)
mcmcplot(out)

# Specify model in BUGS language
Section7p8p4_code.1 <- nimbleCode( {
    
    # Prior distributions
    p0 ~ dunif(0,1)
    alpha0 <- log(p0/(1-p0))
    alpha1 ~ dnorm(0, 0.01)
    alpha2 ~ dnorm(0, 0.01)
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    psi <- sum(lambda[1:50])/M
    tau ~ dgamma(0.1,0.1) # New parameter, precision of ind. random effects
    sigma <- 1/sqrt(tau)
    
    # log-linear model for abundance: lambda depends on WOODY
    for(s in 1:nsites){
      log(lambda[s])<- beta0 + beta1*X[s,1]
      probs[s]<- lambda[s]/sum(lambda[1:50])
    }
    
    # Model for individual encounter histories
    for(i in 1:M){
      eta[i] ~ dnorm(alpha0, tau)  # Individual random effect
      group[i] ~ dcat(probs[1:50])  # Group == site membership
      z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Observation model: p depends on STRUCT + TIME + ind. heterogeneity 
    for(j in 1:J){
      logit(p[i,j]) <- alpha1 * X[group[i],2] + alpha2*X[group[i],3] + eta[i]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
    }   
    }
  }
)

# Parameters monitored: add sigma
parameters <- c("p0", "alpha0", "alpha1", "alpha2", "beta0", "beta1", "psi", "sigma")

# Initial values: add tau
inits <- function(){
    list (p0 = runif(1),
          alpha1 = runif(1),
          alpha2 = rnorm(1),
          beta0=runif(1),    
          beta1=rnorm(1),
          z= c( rep(1,100), rep(0, 300)),
          tau = 1,
# To avoid nimble warnings, we initialize unobserved groups
          group = c(rep(NA,length(alfl$id)),
                    sample(1:nrow(data$X), size = length(site) - length(alfl$id), replace = TRUE))
          )  }

# MCMC settings (others as before)
ni <- 50000   ;   nb <- 10000   ;   nt <- 2   ;   nt <- 10

# Call nimble from R and summarize marginal posteriors
out <- nimbleMCMC(code = Section7p8p4_code.1, 
                  data = data, 
                  constants = constants,
                  inits = inits,
                  monitors = params,
                  nburnin = nb, 
                  niter = ni,
                  samplesAsCodaMCMC  = TRUE)

mcmcplot(out)
# This does not mix extremely well.

# compare relative efficiencies of nimble and jags:
Section7p8p4_compare <- compareMCMCs(
  modelInfo = list(
    code = Section7p8p4_code.1,
    data = data,
    constants = constants,
    inits = inits()
  ),
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = 3*nb,
  niter = 3*ni # Use a longer run.
)

make_MCMC_comparison_pages(Section7p8p4_compare, modelNames = "Section7p8p4", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section7p8p4.html"))
