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
source("Chapter_9_setup.R")
## 9.3.2 Simulating some time-removal/DS data
## ------------------------------------------------------------------------
## Obtain a data set and harvest the results
library(AHMbook)
set.seed(1235)                 ## so we all create the same data set 
temp <- simHDStr(type="point") ## Simulate point count-removal data set
data <- temp$data              ## harvest data
B <- temp$B                    ## upper limit of counting (maximum count distance)
nsites <- temp$nsites          ## Number of sites
habitat <- temp$habitat        ## habitat covariate
K <- temp$K                    ## Number of removal periods

## Create the observed encounter frequencies per site (include the zeros! )
data <- data[!is.na(data[,2]),]   ## Sites where detections did occur
n <- rep(0,nsites)                ## The full site vector
names(n) <- 1:nsites
n[names(table(data[,1]))] <- table(data[,1])  ## Put in the counts
site <- data[,1]
nobs <- nrow(data) 

## Create the distance class data
nD <- 10             ## Number of distance classes 
delta <- B/nD        ## bin size or width
mdpts <- seq(delta/2,B,delta) ## midpoint distance of bins up to max distance
dclass <- data[,"d"] ## distance class for each observation
dclass <- dclass%/%delta  +1
tint <- data[,"aux"]

## Bundle data and summarize
str( win.data <- list(n=n, 
                    site=site, 
                    dclass=as.numeric(dclass),
                    nsites=nsites, 
                    nobs=nobs, 
                    delta=delta, 
                    nD=nD,
                    mdpts=mdpts,
                    B=B, 
                    K=K, 
                    tint=tint, 
                    habitat=habitat) )

library(nimble)
Section9p3p2_code <- nimbleCode({
    ## Prior distributions for basic parameters
    ## Intercepts
    beta.a0 ~ dnorm(0,0.01)    ## intercept for availability
    alpha0 ~ dnorm(0, 0.01)    ## intercept for sigma
    alpha1 ~ dnorm(0,0.01)     ## slope on sigma covariate
    ## Coefficients
    # beta.a1 ~ dnorm(0,0.01)  ## slope for availability covariate, uncomment if including covariates on availability
    beta0 ~ dnorm(0,0.01)      ## intercept for lambda
    beta1 ~dnorm(0,0.01)       ## slope for lambda covariate
    
    for(s in 1:nsites){
      # Add covariates to scale parameter DISTANCE (perceptibility)
      log(sigma[s]) <- alpha0 + alpha1*habitat[s] 
      # Add covariates for availability here TIME-REMOVAL (availability)
      p.a[s] <- exp(beta.a0) / (1+exp(beta.a0)) 
      # Optional covariates on availability
      # exp(beta.a0 + beta.a1*date[s])/(1+exp(beta.a0+beta.a1*date[s]))
      # Distance sampling detection probability model
      for(b in 1:nD){
        log(g[b,s]) <- -mdpts[b]*mdpts[b]/(2*sigma[s]*sigma[s])  # Half-normal  
        f[b,s] <- ( 2*mdpts[b]*delta )/(B*B) # Radial density function
        pi.pd[b,s] <- g[b,s]*f[b,s]  #  Product Pr(detect)*Pr(distribution)
        pi.pd.c[b,s] <- pi.pd[b,s]/pdet[s]  # Conditional probabilities
      }
      pdet[s] <- sum(pi.pd[1:nD,s])  # Probability of detection at all 
      
      # Time-removal probabilities
      for (k in 1:K){
        pi.pa[k,s] <- p.a[s] * pow(1-p.a[s], (k-1))  
        pi.pa.c[k,s] <- pi.pa[k,s]/phi[s] # Conditional probabilities of availability
      }
      phi[s] <- sum(pi.pa[1:K,s]) # Probability of ever available
    }
    # Conditional observation model for categorical covariates
    for(i in 1:nobs){  
      dclass[i] ~ dcat(pi.pd.c[1:nD,site[i]]) 
      tint[i] ~ dcat(pi.pa.c[1:K,site[i]])
    }
    # Abundance model
    for(s in 1:nsites){ 
      # Binomial model for # of captured individuals
      # n[s] ~ dbin(pmarg[s], M[s]) # Formulation b, see text
      # pmarg[s] <- pdet[s]*phi[s] 
      n[s] ~ dbin(pdet[s], N[s])    # Formulation a, see text
      N[s] ~ dbin(phi[s],M[s])      # Number of available individuals
      M[s] ~ dpois(lambda[s])       # Abundance per survey/site/point
      # Add site-level covariates to lambda
      log(lambda[s]) <- beta0 + beta1*habitat[s] 
    }
    # Derived quantities
    Mtot <- sum(M[1:nsites])  # Total population size
    Ntot <- sum(N[1:nsites])  # Total available population size
    PDETmean <- mean(pdet[1:nsites]) # Mean perceptibility across sites
    PHImean <- mean(phi[1:nsites]) # Mean availability across sites
  }
)

## Create initial values (including for M and N) and list parameters to save
Mst <- Nst <- n + 1
inits <- function(){list(M=Mst, 
                         N=Nst, 
                         alpha0=1, 
                         beta0=runif(1,-1,1), 
                         # beta.a1=runif(1,-1,1), ## uncomment if using covariates for availability
                         beta1=runif(1,-1,1), 
                         alpha1=runif(1,-1,1),
                         beta.a0=runif(1,-1,1))}

params <- c("beta.a0", 
            # "beta.a1",   ## uncomment if using covariates for availability
            "alpha0", "alpha1", "beta0", "beta1", "PDETmean", "PHImean", "Mtot", "Ntot")

# MCMC settings
ni <- 50000   ;   nb <- 10000   ;   nt <- 4   ;   nc <- 3

# Run NIMBLE, check convergence and summarize the results
out2a <- nimbleMCMC(code = Section9p3p2_code,
                    constants = win.data,
                    inits = inits,
                    monitors = params,
                    nburnin = nb,
                    niter = ni,
                    samplesAsCodaMCMC = TRUE)

library(mcmcplots)
colnames(out2a)
mcmcplot(out2a)


Section9p3p2_compare <- compareMCMCs(
  modelInfo = list(
    code = Section9p3p2_code,
    data = win.data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = 2*nb,
  niter = 8*ni,
  thin = nt
)

make_MCMC_comparison_pages(Section9p3p2_compare, modelNames = "Section9p3p2", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section9p3p2.html"))
