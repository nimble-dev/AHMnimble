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
# 9.4.1/9.4.2 Simulating MRDS data and analysis in BUGS
# ------------------------------------------------------------------------
# Simulate a double-observer sampling data set
library(AHMbook)
set.seed(1235)
temp <- simHDStr(type="point", method="double") # simulate double observer point count data set
data <- temp$data         # harvest data
B <- temp$B               # upper limit of counting (maximum count distance)
nsites <-temp$nsites      # number of sites
habitat <-temp$habitat    # habitat covariate


# Processing of the data: pad the count vector with 0s etc.
data <- data[!is.na(data[,2]),]
n <- rep(0,nsites)
names(n) <- 1:nsites
n[names(table(data[,1]))] <- table(data[,1])
site <- data[,1]
dclass <- data[,"d"]      # categorical distance class for each observation
aux <- data[,"aux"]       # the auxiliary variable is capture history

# Create the categorical distance variable, use 10 classes here.
nD <- 10
delta <- B/nD # bin width
mdpts <-seq(delta/2,B,delta) # midpoint of bins up to max distance
nobs <- nrow(data) 
dclass <- dclass%/%delta  +1

# Bundle data and look at overview of data
str( win.data <-list(n=n,
                     site=site, 
                     dclass=as.numeric(dclass), 
                     nsites=nsites, 
                     nobs=nobs, 
                     delta=delta, 
                     nD=nD, 
                     mdpts=mdpts, 
                     B=B, 
                     aux=aux, 
                     habitat=habitat) )

library(nimble)
## Define model in BUGS langauge
# "do_model.txt" in book code
Section9p4p2_code <- nimbleCode({
    
    ## Priors for fixed detection parameters
    ## 2 observer detection probability parameters
    logitp1 ~ dnorm(0, 0.01)
    logitp2 ~ dnorm(0, 0.01)
    ## Intercepts
    alpha0 ~ dnorm(0, 0.01)    #3 intercept for sigma
    alpha1 ~ dnorm(0, 0.01)    #3 slope on sigma covariate
    ## Coefficients
    beta0 ~ dnorm(0,0.01)      #3 intercept for lambda
    beta1 ~ dnorm(0,0.01)      #3 slope for lambda covariate
    
    ## Detection scale parameter model
    for(s in 1:nsites){ 
      ## Covariates on scale parameter (perceptibility)
      log(sigma[s]) <- alpha0 + alpha1*habitat[s] 
      ## Double observer cell probabilities here if there are covariates
      logit(pobs[1,s]) <- logitp1 ## + covariates
      logit(pobs[2,s]) <- logitp2 ## + covariates
      
      ## Distance sampling model and cell probabilities
      for(b in 1:nD){
        log(g[b,s]) <- -mdpts[b]*mdpts[b]/(2*sigma[s]*sigma[s])  ## half-normal 
        f[b,s] <- ( 2*mdpts[b]*delta )/(B*B) ## Scaled radial density function
        pi.pd[b,s] <- g[b,s]*f[b,s]          ##  Product Pr(detect)*Pr(distribution)
        pi.pd.c[b,s] <- pi.pd[b,s]/pdet[s]   ## Conditional cell probabilities
      }
      pdet[s] <- sum(pi.pd[1:nD,s])              ## Marginal probability of detection
      
      ## Double observer cell probabilities and conditional probabilities
      doprobs[1,s] <- pobs[1,s]*(1-pobs[2,s]) 
      doprobs.condl[1,s] <- doprobs[1,s]/sum(doprobs[1:3,s])
      doprobs[2,s] <- (1-pobs[1,s])*pobs[2,s]
      doprobs.condl[2,s] <- doprobs[2,s]/sum(doprobs[1:3,s])
      doprobs[3,s] <- pobs[1,s]*pobs[2,s] 
      doprobs.condl[3,s] <- doprobs[3,s]/sum(doprobs[1:3,s])
      pavail[s] <- sum(doprobs[1:3,s])  ## probability of availability AT ALL
    }
    
    ## Observation model for two categorical covariates
    for(i in 1:nobs){  
      dclass[i] ~ dcat(pi.pd.c[1:nD,site[i]]) 
      aux[i] ~ dcat(doprobs.condl[1:3,site[i]])
    }
    
    ## Abundance model
    for(s in 1:nsites){ 
      ## Binomial model for # of captured individuals
      n[s] ~ dbin(pdet[s], N[s]) 
      N[s] ~ dbin(pavail[s], M[s])   ## binomial availability model
      ## Abundance model
      M[s] ~ dpois(lambda[s])        ## predicted abundance per survey/site/point
      ## Add site-level covariates to lambda
      log(lambda[s])<- beta0 + beta1*habitat[s] 
      }
    ## Derived parameters
    Mtot <- sum(M[1:nsites])
    Ntot <- sum(N[1:nsites])
    logit(p1) <- logitp1
    logit(p2) <- logitp2
    sigma0 <- exp(alpha0)            ## Baseline sigma
  }
)

## Inits function
Nst <- n + 1         ## inits for N
inits <- function(){list(M=Nst+1, 
                         N=Nst, 
                         alpha0=runif(1,1,2),
                         beta0=runif(1,-1,1), 
                         beta1=runif(1,-1,1), 
                         alpha1=runif(1,-1,1),
                         logitp1=0, 
                         logitp2=0)} 

# Parameters to monitor
params <- c("alpha0", "alpha1", "beta0", "beta1", "Ntot", "Mtot", "logitp1", 
   "logitp2", "p1", "p2", "sigma0")

# MCMC settings
ni <- 50000   ;   nb <- 10000   ;   nt <- 4   ;   nc <- 3

# Run NIMBLE, check convergence and summarize the results
out3 <- nimbleMCMC(
  code = Section9p4p2_code,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

## look at posterior information
library(mcmcplots)
colnames(out3)
mcmcplot(out3)

## Compare MCMC efficiency of nimble to jags
Section9p4p2_compare <- compareMCMCs(
  modelInfo = list(
    code = Section9p4p2_code,
    data = win.data,
    inits = inits()
  ),
  monitors = params,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = nb, 
  niter = 2*ni
)

make_MCMC_comparison_pages(Section9p4p2_compare, modelNames = "Section9p4p2", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section9p4p2.html"))
