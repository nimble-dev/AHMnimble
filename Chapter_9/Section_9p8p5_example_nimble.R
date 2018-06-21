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
## This example builds slowly in nimble.  It may go faster
## with a deeper re-write of the model.
source("Chapter_9_setup.R")
# 9.8.5 Hierarchical spatial distance sampling
# ------------------------------------------------------------------------
library(AHMbook)
#sim.spatialHDS(lam0 = 4, sigma = 1.5, B = 3, nsites = 100)
# lam0 = expected population size per site
# nsites = number of point count locations
# B = count radius. Function simulates coordinates of individuals on a square
#       [0,2*B] x[0,2*B], with a count location on the point (B,B)
# sigma = scale of half-normal detection function

library(raster)

# Simulate a data set and harvest the output
set.seed(1234)
str(tmp <-sim.spatialHDS(lam0 = 3, sigma = 1.5, B = 3, nsites = 100, show.plot = 0))

# Process the simulated data set
data <- tmp$data
# To make it a �real� data set:
data <- data[!is.na(data[,2]),] # Get rid of the 0 sites 
data <- data[data[,"y"]==1,]    # Only keep detected individuals

# Now zero-pad the observed counts
nsites <- tmp$nsites
nobs <- rep(0, nsites)
names(nobs) <- 1:nsites
nobs[names(table(data[,1]))] <- table(data[,1])

# Extract data elements that we need
site <- data[,"site"]
s <- data[,c("u1","u2")]
B <- tmp$B
Habitat <- (tmp$Habitat) # Raster values
Habgrid <- tmp$grid   # Raster coordinates
nind <- nrow(data)
G <- nrow(Habgrid)

# We have to convert observed locations to pixels
pixel <- rep(NA,nrow(data))
D <- e2dist(s[1:nind,], Habgrid)
for(i in 1:nind){
  pixel[i] <- (1:ncol(D))[D[i,]==min(D[i,])]
}

# Do data augmentation of data from each site �S-fold� DA. Three objects need
# to have DA applied to them: Ymat, umat, pixmat
Msite <- 2*max(nobs)  # Perhaps use a larger value 
Ymat <- matrix(0,nrow=Msite,ncol=nsites)
umat <- array(NA, dim=c(Msite,nsites,2))
pixmat <- matrix(NA,nrow=Msite,ncol=nsites)
for(i in 1:nsites){
  if(nobs[i]==0) next
  Ymat[1:nobs[i],i]<- data[data[,1]==i,"y"]
  umat[1:nobs[i],i,1:2]<- data[data[,1]==i,c("u1","u2")]
  pixmat[1:nobs[i],i]<- pixel[data[,1]==i]
}

# Bundle the data for BUGS
str(data <- list (y=Ymat, pixel=pixmat, Habitat=Habitat, Habgrid=Habgrid, G = G, 
                  nsites=nsites, M = Msite, B = B))

library(nimble)
# Write out the BUGS model file
# "spatialHDS.txt" in book code
Section9p8p5_code <- nimbleCode({ 
    
    # Prior distributions
    sigma ~ dunif(0,10)
    beta1 ~ dnorm(0,0.01)
    beta0 ~ dnorm(0,0.01)
    lam0 <- exp(beta0)*G           # Baseline lambda in terms of E(N) per sample unit
    
    # For each site, construct the DA parameter as a function of lambda
    for(s in 1:nsites){
      lamT[s] <- sum(lambda[1:G,s])   # total abundance at a site
      psi[s] <- lamT[s]/M
      for(g in 1:G){             # g is the pixel index, there are G total pixels
        lambda[g,s] <- exp(beta0 + beta1*Habitat[g,s])
        probs[g,s] <- lambda[g,s]/sum(lambda[1:G,s])
      }
    }
    
    # DA variables and spatial location variables:
    for(s in 1:nsites){
      for(i in 1:M){
        z[i,s] ~ dbern(psi[s])
        pixel[i,s] ~ dcat(probs[1:G,s])
        u[i,s,1:2] <- Habgrid[pixel[i,s],]   # location = derived quantity  
        # distance = derived quantity
        d[i,s] <- pow(   pow( u[i,s,1]-B,2) + pow(u[i,s,2]-B,2), 0.5)
        p[i,s] <- exp(-d[i,s]*d[i,s]/(2*sigma*sigma))  # Half-normal model
        mu[i,s] <- p[i,s]*z[i,s]
        y[i,s] ~ dbern(mu[i,s])    # Observation model
      }
      # Derived parameters
      N[s]<- sum(z[1:M,s])            # Site specific abundance
    }
    Ntotal <- sum(N[1:nsites])             # Total across all sites
    D <- Ntotal/(9*nsites)         # Density: point area = 9 ha 
  }
)


# Inits and parameters saved
zst <- Ymat
inits <- function(){
    # Add inits for pixel to avoid nimble warnings
    pixelInit <- data$pixel
    pixelIsData <- !is.na(pixelInit)
    pixelInit[pixelIsData] <- NA
    pixelInit[!pixelIsData] <- sample(1:G, size = sum(!pixelIsData), replace = TRUE)
    list (sigma=1.5,
          z = zst,
          beta0 = -5,
          beta1=1,
          pixel = pixelInit) }


params <- c("sigma", "Ntotal", "beta1", "beta0", "D", "lam0")

# The stochastic indexing of Habgrid requires that we separate constants
# from data
constants <- data[c('B', 'M', 'nsites', 'G')]
data <- data[c('y', 'Habitat', 'Habgrid', 'pixel')]
data$Habgrid <- as.matrix(data$Habgrid)


# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nthin <- 2   ;   nc <- 3

# Call nimble, check convergence and summarize the results
out3 <- nimbleMCMC(
  code = Section9p8p5_code,
  constants = constants,
  data = data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

library(mcmcplots)
library(nimble)
library(coda)
colnames(out3)
mcmcplot(out3)

initsValues <- inits()

Section9p8p5_compare <- compareMCMCs(
  modelInfo = list(
    code = Section9p8p5_code,
    data = data,
    constants = constants,
    inits = initsValues
  ),
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  monitors = params,
  burnin = 2*nb,
  niter = 2*ni
)

make_MCMC_comparison_pages(Section9p8p5_compare, modelNames = "Section9p8p5", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section9p8p5.html"))

