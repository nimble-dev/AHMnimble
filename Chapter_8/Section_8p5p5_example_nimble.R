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
# 8.5.5 Analysis of the ISSJ data
# ------------------------------------------------------------------------
# Load the ISSJ data
library(unmarked)
data(issj)

# Prepare some data
nD <- 3                          # Number of intervals
delta <- 100                     # Interval width
B <- 300                         # Upper bound (max. distance)
midpt <- c(50, 150, 250)         # mid points

# Convert vector frequencies to individual distance class
H <- as.matrix(issj[,1:3])
nsites <- nrow(H)
ncap <- apply(H, 1, sum)         # Number of individuals detected per site
dclass <- rep(col(H), H)         # Distance class of each individual
nind <- length(dclass)           # Number of individuals detected
elevation <- as.vector(scale(issj[,c("elevation")])) # Prepare covariates
forest <- as.vector(scale(issj[,"forest"]))
chaparral <- as.vector(scale(issj[,"chaparral"]))


# Bundle and summarize data set
str( win.data <- list(nsites=nsites, 
                      nind=nind, 
                      B=B, 
                      nD=nD, 
                      midpt=midpt,
                      delta=delta, 
                      ncap=ncap, 
                      chaparral=chaparral, 
                      elevation=elevation, 
                      dclass=dclass) )

require(nimble)

# BUGS model specification
# "model5.txt" in the book code
Section8p5p5_code <- nimbleCode({ 
    # Priors
    sigma ~ dunif(0,1000)
    beta0 ~ dunif(-10,10)
    beta1 ~ dunif(-10,10)
    beta2 ~ dunif(-10,10)
    beta3 ~ dunif(-10,10)
    sigma.site ~ dunif(0,10)
    tau <- 1/(sigma.site*sigma.site)
    
    # Specify hierarchical model
    for(i in 1:nind){
      dclass[i] ~ dcat(fc[1:nD]) # Part 1 of HM
    }
    
    # construct cell probabilities for nG cells
    for(g in 1:nD){                # midpt = mid-point of each cell
      log(p[g]) <- -midpt[g] * midpt[g] / (2 * sigma * sigma)
      pi[g] <- ((2 * midpt[g]) / (B * B)) * delta # prob. per interval
      f[g] <- p[g] * pi[g]
      fc[g] <- f[g] / pcap
    }
    pcap <- sum(f[1:nD])               # Pr(capture): sum of rectangular areas
    
    for(s in 1:nsites){
    ncap[s] ~ dbin(pcap, N[s])   # Part 2 of HM
    N[s] ~ dpois(lambda[s])      # Part 3 of HM
    log(lambda[s]) <- beta0 + beta1*elevation[s] + beta2*chaparral[s] + beta3*chaparral[s]*chaparral[s] + site.eff[s] # linear model for abundance
    site.eff[s] ~ dnorm(0, tau)   # Site log-normal 'residuals'
    }
    # Derived params
    Ntotal <- sum(N[1:nsites])
    area<- nsites*3.141*300*300/10000   # Total area sampled, ha
    D<- Ntotal/area
  }
)


# Inits
Nst <- ncap + 1
inits <- function(){list (sigma = runif(1, 30, 100), 
                          beta0 = 0, 
                          beta1 = 0, 
                          beta2 = 0, 
                          beta3 = 0, 
                          N = Nst, 
                          sigma.site = 0.2)}

# Params to save
params <- c("sigma", "beta0", "beta1", "beta2", "beta3", "sigma.site", "Ntotal","D")

# MCMC settings
ni <- 52000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3

# Run nimble from R and summarize posteriors
out5 <- nimbleMCMC(code = Section8p5p5_code, 
                   constants = win.data, 
                   inits = inits, 
                   monitors = params,
                   nburnin = nb, 
                   niter = ni,
                   samplesAsCodaMCMC = TRUE) 

## examine posterior output
library(coda)
library(mcmcplots)
colnames(out5)
mcmcplot(out5)
## all beta parameters exhibit slow mixing

Section8p5p5_compare <- compareMCMCs(
  modelInfo = list(
    code = Section8p5p5_code,
    data = win.data,
    inits = inits()
  ),
  MCMCs <- c('nimble', 'jags'),
  monitors = params,
  summary = FALSE,
  burnin = nb,
  niter = 2*ni # Larger sample size to improve accuracy of ESS estimates
)

make_MCMC_comparison_pages(Section8p5p5_compare, modelNames = "Section8p5p5", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section8p5p5.html"))

