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
# 8.5.4 Point transect HDS using the conditional multinomial formulation
# ------------------------------------------------------------------------
# Simulate a data set using our simHDS function
library(AHMbook)
set.seed(1234)
tmp <- simHDS(type="point", discard0=FALSE)
attach(tmp)

# Prepare data
# Number of individuals detected per site
ncap <- table(data[,1])            # ncap = 1 if no individuals captured
sites0 <- data[is.na(data[,2]),][,1] # sites where nothing was seen
ncap[as.character(sites0)] <- 0    # Fill in 0 for sites with no detections
ncap <- as.vector(ncap)            # Number of individuals detected per site

# Other data
site <- data[!is.na(data[,2]),1]   # Site ID of each observation
delta <- 0.1                       # Distance bin width for rect. approx.
midpt <- seq(delta/2, B, delta)    # Make mid-points and chop up data
dclass <- data[,5] %/% delta + 1   # Convert distance to distance category
nD <- length(midpt)                # Number of distance intervals
dclass <- dclass[!is.na(data[,2])] # Observed categorical observations
nind <- length(dclass)             # Total number of individuals detected

# Bundle and summarize data set
str( win.data <- list(nsites=nsites, 
                      nind=nind, 
                      B=B, 
                      nD=nD, 
                      midpt=midpt,
                      delta=delta, 
                      ncap=ncap, 
                      habitat=habitat, 
                      wind=wind, 
                      dclass=dclass,
                      site=site) )

require(nimble)
# BUGS model specification for point transect data
# "model4.txt" in book code
Section8p5p4_code <- nimbleCode({ 
    # Priors
    alpha0 ~ dunif(-10,10)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-10,10)
    beta1 ~ dunif(-10,10)
    
    for(i in 1:nind){
      dclass[i] ~ dcat(fc[site[i],1:nD]) # Part 1 of HM
    }
    for(s in 1:nsites){
      # Construct cell probabilities for nD distance bands
      for(g in 1:nD){                # midpt = mid-point of each band
        log(p[s,g]) <- -midpt[g] * midpt[g] / (2 * sigma[s] * sigma[s])
        pi[s,g] <- ((2 * midpt[g] ) / (B * B)) * delta # prob. per interval
        f[s,g] <- p[s,g] * pi[s,g]
        fc[s,g] <- f[s,g] / pcap[s]
      }
      pcap[s] <- sum(f[s,1:nD])           # Pr(capture): sum of rectangular areas
      
      ncap[s] ~ dbin(pcap[s], N[s])   # Part 2 of HM
      N[s] ~ dpois(lambda[s])         # Part 3 of HM
      log(lambda[s]) <- beta0 + beta1 * habitat[s] # linear model abundance
      log(sigma[s]) <- alpha0 + alpha1*wind[s]     # linear model detection
    }
    
    # Derived parameters
    Ntotal <- sum(N[1:nsites])
    area <- nsites*3.141*B*B
    D <- Ntotal/area
  }
)

# Inits
Nst <- ncap + 1
inits <- function(){list(alpha0=0, 
                         alpha1=0, 
                         beta0=0, 
                         beta1=0, 
                         N=Nst)}

# Params to save
params <- c("alpha0", "alpha1", "beta0", "beta1", "Ntotal","D")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 1   ;   nc <- 3

# Run nimble from R
out4 <- nimbleMCMC(code = Section8p5p4_code, 
                   constants = win.data, 
                   inits = inits, 
                   monitors = params,
                   nburnin = nb, 
                   niter = ni,
                   samplesAsCodaMCMC = TRUE)

## examine posterior informaton
library(coda)
library(mcmcplots)
colnames(out4)
mcmcplot(out4)

Section8p5p4_compare = compareMCMCs(
  modelInfo = list(
    code = Section8p5p4_code,
    data = win.data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = 3*nb,
  niter = 3*ni
)

make_MCMC_comparison_pages(Section8p5p4_compare, modelNames = "Section8p5p4", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section8p5p4.html"))
