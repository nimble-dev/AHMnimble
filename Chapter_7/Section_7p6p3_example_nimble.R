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
source("Section_7p6_setup.R")
# 7.6.3 Poisson formulation of the multinomial mixture model
# ------------------------------------------------------------------------
require(nimble)
# Specify model in BUGS language
# This code is called "modelP.txt" in original book code.
Section7p6p3_code <- nimbleCode( {
    
    # Prior distributions
    p0 ~ dunif(0,1)
    alpha0 <- logit(p0)
    alpha1 ~ dnorm(0, 0.01)
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    beta2 ~ dnorm(0, 0.01)
    beta3 ~ dnorm(0, 0.01)
    
    for(i in 1:M){
      # logit-linear model for detection: understory cover effect
      logit(p[i]) <- alpha0 + alpha1 * X[i,1]
      # log-linear model for abundance: UFC + TRBA + UFC:TRBA
      log(lambda[i])<- beta0 + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,2]*X[i,1]
    
    # Poisson parameter = multinomial cellprobs x expected abundance
    pi[i,1] <- p[i] * lambda[i]
    pi[i,2] <- p[i] * (1-p[i]) * lambda[i]
    pi[i,3] <- p[i] * (1-p[i]) * (1-p[i]) * lambda[i]
    pi[i,4] <- p[i] * (1-p[i]) * (1-p[i]) * (1-p[i]) * lambda[i]
    
    for(j in 1:4){
      y[i,j] ~ dpois(pi[i,j])
    }
      # Generate predictions of N[i]
      N[i] ~ dpois(lambda[i])
    }
  }
)
inits <- function(){
    list (p0 = runif(1),
          alpha1=runif(1),
          beta0=runif(1),
          beta1=runif(1),
          beta2=runif(1),
          beta3=runif(1))
}

# Define parameters to save and MCMC settings
params <- c("p0", "alpha1", "beta0", "beta1", "beta2", "beta3", "N")
nc <- 3   ;   ni <- 6000   ;   nb <- 1000   ;   nt <- 1

# Call nimble from R and summarize marginal posteriors 
out <- nimbleMCMC(code = Section7p6p3_code, 
                  constants = data, 
                  inits = inits,
                  monitors = params,
                  nburnin = nb, 
                  niter = ni,
                  samplesAsCodaMCMC = TRUE)

## look at posterior information
require(coda)
require(mcmcplots)
colnames(out)
params_of_interest <- c("alpha1", "beta0", "beta1", "beta2", "beta3", "p0")
mcmcplot(out[, params_of_interest])

# Specify model in BUGS language
# This code is also called "modelP.txt" in original book code.
Section7p6p3_code.1 <- nimbleCode( {
    
    # Prior distributions
    p0 ~ dunif(0,1)
    alpha0 <- logit(p0)
    alpha1 ~ dnorm(0, 0.01)
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    beta2 ~ dnorm(0, 0.01)
    beta3 ~ dnorm(0, 0.01)
    
    tau ~ dgamma(0.1,0.1)  # Excess-Poisson variation (precision)
    sigma <- sqrt(1 / tau)
    
    for(i in 1:M){
      # logit-linear model for detection: understory cover effect
      logit(p[i]) <- alpha0 + alpha1 * X[i,1]
      # Normal random effects
      eta[i] ~ dnorm(0, tau)  # 'residuals' for extra-Poisson noise
      # log-linear model for abundance: UFC + TRBA + UFC:TRBA + eta
      log(lambda[i])<- beta0 + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,2]*X[i,1] + eta[i]     # note 'extra-residual' for overdispersion
    
    # Poisson parameter = multinomial cellprobs x expected abundance
    pi[i,1] <- p[i] * lambda[i]
    pi[i,2] <- p[i] * (1-p[i]) * lambda[i]
    pi[i,3] <- p[i] * (1-p[i]) * (1-p[i]) * lambda[i]
    pi[i,4] <- p[i] * (1-p[i]) * (1-p[i]) * (1-p[i]) * lambda[i]
    
    for(j in 1:4){
      y[i,j] ~ dpois(pi[i,j])
    }
    # Generate predictions of N[i]
    N[i] ~ dpois(lambda[i])
    }
  }
)

# Inits
inits <- function(){
  list (p0 = runif(1), 
        alpha1=runif(1), 
        beta0=runif(1), 
        beta1=runif(1), 
         beta2=runif(1), 
        beta3=runif(1), 
        tau = 1)}

# Define parameters to save and MCMC settings
params <- c("p0", "alpha1", "beta0", "beta1", "beta2", "beta3", "sigma")
nc <- 3   ;   ni <- 32000   ;   nb <- 2000   ;   nt <- 1

# Call nimble from R and summarize marginal posteriors
out <- nimbleMCMC(code = Section7p6p3_code.1, 
                  constants = data, 
                  inits = inits, 
                  monitors = params,
                  nburnin = nb, 
                  niter = ni,
                  samplesAsCodaMCMC = TRUE)

## look at posterior information
colnames(out)
params_of_interest <- c("p0", "alpha1", "beta0", "beta1", "beta2", "beta3")
mcmcplot(out[, params_of_interest])
## samples appear to be well-mixed based on visual inspection of traceplot and autocorrelation

## compare mcmc efficiencies
Section7p6p3_compare <- compareMCMCs(
  modelInfo = list(
    code = Section7p6p3_code.1,
    data = data,
    inits = inits(),
    monitors = params),
  MCMCs = c('jags', 'nimble'),
  summary = FALSE,
  burnin = nb,
  niter = ni
)

make_MCMC_comparison_pages(Section7p6p3_compare, modelNames = "Section7p6p3", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section7p6p3.html"))


