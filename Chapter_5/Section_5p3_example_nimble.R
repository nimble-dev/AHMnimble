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
## NEEDS WORK
source("Chapter_5_setup.R")

require(nimble)
# Code for the file "multiple_linear_regression_model.txt"
# from the original book code.
section5p3_code <- nimbleCode({

                                        # Priors
    alpha0 ~ dnorm(0, 1.0E-06)           # Prior for intercept
    alpha1 ~ dnorm(0, 1.0E-06)           # Prior for slope of elev
    alpha2 ~ dnorm(0, 1.0E-06)           # Prior for slope of forest
    alpha3 ~ dnorm(0, 1.0E-06)           # Prior for slope of interaction
    tau <- pow(sd, -2)                   # Precision tau = 1/(sd^2)
    sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale

                                        # Likelihood
    for (i in 1:M){
        Cmean[i] ~ dnorm(mu[i], tau)      # dispersion tau is precision (1/variance)
        mu[i] <- alpha0 + alpha1*elev[i] + alpha2*forest[i] + alpha3*elev[i]*forest[i]
    }

                                        # Derived quantities
    for (i in 1:M){
        resi[i] <- Cmean[i] - mu[i]
    }
})

## BUGS and JAGS for "inits" (initial values) expect a *function*
## NIMBLE takes a list of values 
## So, for now you'll need to call the function.

inits <- function() list(alpha0 = rnorm(1, 0, 10), alpha1 = rnorm(1, 0, 10),    ## Define inits function (not in original code)
                         alpha2 = rnorm(1, 0, 10), alpha3 = rnorm(1, 0, 10))

initsValues <- inits()

# Parameters monitored (i.e., for which estimates are saved)
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "resi")

## Guesses at what these are: nb is burnin; ni is number of iterations to save; nc is number of chains (nimble only supports 1!); nt is thinning interval.
# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3
# ni <- 10   ;   nt <- 1   ;   nb <- 0   ;  nc <- 8 # not run

## Get MCMC running in nimble.  For illustration we will write out all the steps here rather than use the nimbleMCMC function.

## 1. make model

section5p3_model <- nimbleModel(section5p3_code,
                                constants = win.data,
                                inits = initsValues) ## sometimes we need dimensions
## for big models, make it faster with check = FALSE and/or calculate = FALSE

## 2. configure MCMC
section5p3_mcmcConf <- configureMCMC(section5p3_model, monitors = params)
section5p3_mcmcConf$printSamplers()

## 3. build MCMC
section5p3_mcmc <- buildMCMC(section5p3_mcmcConf)

## 4. compile model and MCMC (can do this in one step or two)
Csection5p3_model <- compileNimble(section5p3_model)
Csection5p3_mcmc <- compileNimble(section5p3_mcmc, project = section5p3_model)

## 5a. Run MCMC for burnin
Csection5p3_samples <- runMCMC(Csection5p3_mcmc, niter = ni, nburnin = nb)
Csection5p3_mcmc$run(niter = nb)

## 5b. Run MCMC to save samples
Csection5p3_mcmc$run(niter = ni)

## 6. extract samples
Csection5p3_samples <- as.matrix(Csection5p3_mcmc$mvSamples)

require(coda)
## 7. convert samples to an mcmc object for the coda package
out1N <- as.mcmc(Csection5p3_samples)

## AHM book code uses jagsUI::traceplot.  jagsUI does not provide an easy function
## for converting other MCMC outputs to one of their jagsUI objects, so we cannot
## easily use jagsUI::traceplot.  Instead we use package mcmcplots to look at the same
## parameters as chosen in AHM.

library(mcmcplots)
mcmcplot(out1N[, c('alpha1','alpha2', Csection5p3_model$expandNodeNames('resi[c(1,3, 5:6)]'))])

## plot and summarize output
plot(out1N)
summary(out1N)

## Next we can run NIMBLE and JAGS on the same model to compare efficiency.
## We define MCMC efficiency as (effective sample size) / (computation time).
## We pay particular attention to the minimum (over parameters) MCMC efficiency,
## as that reflects the worst-mixing parameter.
Section5p3_compare <- compareMCMCs(
  modelInfo = list(code = section5p3_code,
                  data = win.data, 
                  inits = initsValues),
  MCMCs = c('nimble', 'jags', 'nimble_slice'),
  summary = TRUE,
  burnin = nb,
  niter = ni)
    
make_MCMC_comparison_pages(Section5p3_compare,
                           modelNames = "Section5p3",
                           dir = outputDirectory)

## This should open the results in a browser.
browseURL(Sys.glob(file.path(outputDirectory, 'Section5p3.html')))
## If it fails, simply use a browser to open the Section5p3.html file
