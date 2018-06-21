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
source("Chapter8_dataSimulation.R")
# 8.3.4 Bayesian analysis of point transect data
# ------------------------------------------------------------------------
### Version 1: Point count data in BUGS (conditional likelihood)
# Simulate a data set and harvest the output
set.seed(1234)
tmp <- sim.pdata(N=200, sigma=1, keep.all=FALSE, B=3)
attach(tmp)

# Chop the data into bins
delta <- 0.1           # width of distance bins for approximation
xg <- seq(0, B, delta) # Make the mid points and chop up the data
midpt <- xg[-1] - delta/2

# Convert distances to categorical distances (which bin?)
dclass <- d %/% delta + 1
nD <- length(midpt) # how many intervals
nind <- length(dclass)

# Bundle and summarize data set
str(win.data <- list(midpt=midpt, 
                     delta=delta, 
                     B=B, 
                     nind=nind, 
                     nD=nD, 
                     dclass=dclass) )


require(nimble)
# BUGS model specification, conditional version
# "model3.txt" in original book code
Section8p3p4_code_model3 <- nimbleCode({
    
    # Prior for single parameter
    sigma ~ dunif(0, 10)
    
    # Construct cell probabilities for nG cells (rectangle approximation)
    for(g in 1:nD){    # midpt[g] = midpoint of each distance band
      log(p[g]) <- -midpt[g] * midpt[g] / (2*sigma*sigma)
      pi[g] <- (( 2 * midpt[g] ) / (B*B)) * delta
      f[g] <- p[g] * pi[g]
      fc[g] <- f[g] / pcap
    }
    pcap <- sum(f[1:nD]) # capture prob. is the sum of all rectangular areas
    
    # Categorical observation model
    for(i in 1:nind){
      dclass[i] ~ dcat(fc[1:nD])
    }
    # Derived quantity: population size
    N <- nind / pcap 
    D<- N/(3.141*B*B)
  }
)

# Inits function
inits <- function(){list (sigma=runif(1, 1, 10)) }

# Params to save
params <- c("sigma", "N","D")

# MCMC settings
ni <- 62000   ;   nb <- 2000   ;   nt   <-   2   ;   nc <- 3

# Run BUGS and summarize posteriors
out3 <- nimbleMCMC(code = Section8p3p4_code_model3, 
                   constants = win.data, 
                   inits = inits, 
                   monitors = params,
                   nburnin = nb, 
                   niter = ni,
                   samplesAsCodaMCMC = TRUE)

## look at posterior information
require(coda)
require(mcmcplots)
mcmcplot(out3)
## mixing looks good based on visual inspection of traceplot and autocorrelation

## compare MCMC efficiencies
Section8p3p4_compare_model3 <- compareMCMCs(
  modelInfo = list(
    code = Section8p3p4_code_model3,
    data = win.data,
    inits = inits()
  ),
  MCMCs = c('nimble', 'jags'),
  monitors = params,
  summary = FALSE,
  burnin = nb,
  niter = ni
)

make_MCMC_comparison_pages(Section8p3p4_compare_model3, modelNames = "Section8p3p4_model3", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section8p3p4_model3.html"))

## Version 2: point count data (full likelihood with data augmentation)
# Do data augmentation (for same simulated data set)
M <- 400
nz <- M - nind
y <- c(rep(1, nind), rep(0, nz))
dclass <- c(dclass, rep(NA, nz))

# Bundle and summarize data set
str( win.data <- list(midpt=midpt, 
                      delta=delta, 
                      B=B, 
                      nind=nind, 
                      nD=nD, 
                      dclass=dclass, 
                      y=y, 
                      nz=nz) )

# BUGS model
# "model4.txt" in book code
Section8p3p4_code_model4 <- nimbleCode({ 
    
    # Priors
    sigma ~ dunif(0, 10)
    psi ~ dunif(0, 1)
    
    # Construct cell probabilities for nG cells (rectangle approximation)
    for(g in 1:nD){           # midpt[g] = midpoint of each distance band
      log(p[g]) <- -midpt[g] * midpt[g] / (2*sigma*sigma)
      pi[g] <- ((2 * midpt[g]) / (B * B)) * delta
      pi.probs[g] <- pi[g] / norm
      f[g] <- p[g] * pi[g]
      fc[g] <- f[g] / pcap   # conditional probabilities
    }
    pcap <- sum(f[1:nD])# capture prob. is the sum of all rectangular areas
    norm <- sum(pi[1:nD])
    
    # Categorical observation model
    for(i in 1:(nind+nz)){
      z[i] ~ dbern(psi)
      dclass[i] ~ dcat(pi.probs[1:nD])
      mu[i] <- p[dclass[i]] * z[i]
      y[i] ~ dbern(mu[i])
    }
    
    # Derived quantity: population size
    N <- sum(z[1:(nind+nz)])
    D<- N/(3.141*B*B)
  }
)

# Inits
inits <- function(){list (sigma=runif(1,1,10), 
                          psi=runif(1),
                          # inits for z and dtclass are added for NIMBLE
                          z = c(
                              rep(1, win.data$nind),
                              sample(c(0,1), win.data$nz, replace = TRUE, prob = c(0.9, 0.1))
                          ),
                          dclass = c(
                              rep(NA, win.data$nind),
                              sample(1:win.data$nD, win.data$nz, replace = TRUE ))
                          )}

# Params to save
params <- c("sigma", "N","D","psi")

# MCMC settings
ni <- 62000   ;   nb <- 2000   ;   nt   <-   2   ;   nc <- 3

# This model will run very slowly in NIMBLE and JAGS
# due to the large dcat distribution.

# Run BUGS and summarize posteriors
out4 <- nimbleMCMC(code = Section8p3p4_code_model4, 
                   constants = win.data, 
                   inits = inits, 
                   monitors = params,
                   nburnin = nb, 
                   niter = ni,
                   samplesAsCodaMCMC = TRUE)

colnames(out4)
mcmcplot(out4)


# Compare efficiencies between nimble and jags:
Section8p3p4_compare_model4 <- compareMCMCs(
  modelInfo = list(
    code = Section8p3p4_code_model4,
    data = win.data,
    inits = inits()
    ),
  monitors = params,  
  MCMCs = c('nimble', 'jags'),
  burnin = nb,
  niter = ni
)

make_MCMC_comparison_pages(Section8p3p4_compare_model4, modelNames = "Section8p3p4_model4", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section8p3p4_model4.html"))
