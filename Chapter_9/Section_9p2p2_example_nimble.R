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
source("Section_9p2p1_setup.R")
# 9.2.2 Analysis in BUGS
# ------------------------------------------------------------------------

set.seed(1234)                 # we all create same data set
temp <- simHDSg(type="line")   # Execute function
data <- temp$data              # harvest data
B <- temp$B                    # Get strip half width
habitat <- temp$habitat        # habitat covariate
nsites <- temp$nsites          # Number of spatial replicates
groupsize <- data[,"gs"] -1    # Input groupsize-1 as data

M <- 400                        # Size of augmented data set is M
nz <- M-nrow(data)              # Number of "pseudo-groups" added
y <- c(data[,2],rep(0,nz))      # Indicator of capture (== 1 for all obs. groups)
nind <- nrow(data)              # Number of observed groups
site <- c(data[,1], rep(NA,nz)) # Site they belong to is unknown 
d <- c(data[,5], rep(NA,nz))    # Their distance data are missing ...
groupsize <- c(groupsize, rep(NA,nz)) # .... as is their size
zst <- y                        # Starting values for data augmentation variable

## Bundle data and produce summary
str(win.data <- list (y=y, B=B, nind=nind, nsites=nsites, d=d, habitat=habitat, 
                       site=site, nz=nz, groupsize=groupsize))

library(nimble)
## Define model in BUGS langauge as nimbleCode object
# "model1.txt" in book code
Section9p2p2_code <- nimbleCode({ 
    
    ## Prior distributions for model parameters
    alpha0 ~ dunif(-10,10)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-10,10)
    beta1 ~ dunif(-10,10)
    lambda.group ~ dgamma(0.1, 0.1)
    
    ## psi is a derived parameter
    psi <- sum(lambda[1:nsites])/(nind+nz)
    
    ## Individual level model: observations and process
    for(i in 1:(nind+nz)){
      z[i] ~ dbern(psi)                   # Data augmentation variables
      d[i] ~ dunif(0, B)                  # Distance is uniformly distributed
      groupsize[i] ~ dpois(lambda.group)  # Group size is Poisson
      
      log(sigma[i]) <- alpha0 +  alpha1*groupsize[i]
      mu[i] <- z[i]*exp(-d[i]*d[i]/(2*sigma[i]*sigma[i])) #p dep on dist class 
      
      ## here using the half normal detection function (Buckland et al. 2001)
      y[i] ~ dbern(mu[i])
      
      site[i] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
      zg[i]<- z[i]*(groupsize[i] + 1)      # Number of individuals in that group
    }
    
    for(s in 1:nsites){
      ## Model for population size of groups
      N[s] ~ dpois(lambda[s])
      log(lambda[s])<- beta0 + beta1*habitat[s]
      site.probs[s]<- lambda[s]/sum(lambda[1:nsites])
    }
    
    # Derived quantities
    G <- sum(z[1:(nind+nz)])        # Total number of groups
    Ntotal <- sum(zg[1:(nind+nz)])  # Total population size (all groups combined)
  }
)

# define MCMC settings, inits function and parameters to save
ni <- 6000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3 ## MCMC settings
inits <- function(){list(alpha0=0,     ## Inits as a function
                         alpha1=0.5, 
                         beta0=0, 
                         beta1=0, 
                         z=zst)}

params <- c("alpha0", "alpha1", "beta0", "beta1", "psi", "Ntotal", "G", 
   "lambda.group")                     ## define params for monitors

# Call NIMBLE, plot posterior distributions
out1 <- nimbleMCMC(code = Section9p2p2_code, 
                   constants = win.data,
                   inits = inits,
                   monitors = params,
                   niter = ni,
                   nburnin = nb,
                   samplesAsCodaMCMC = TRUE)

library(mcmcplots)
colnames(out1)
mcmcplot(out1)
## Mixing is ok but not great.

## Compare nimble efficiency to jags:
Section9p2p2_compare <- compareMCMCs(
  modelInfo = list(
    code = Section9p2p2_code,
    data = win.data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = nb,
  niter = 4*ni  # User larger sample size for more accurate ESS estimates.
)

make_MCMC_comparison_pages(Section9p2p2_compare, modelNames = "Section9p2p2", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section9p2p2.html"))

