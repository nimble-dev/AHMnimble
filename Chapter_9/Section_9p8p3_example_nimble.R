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
# 9.8.3 Modelling spatial covariates
# ------------------------------------------------------------------------
# Simulator function for spatial distance sampling data
library(AHMbook)
sim.spatialDS <- 
  function(N=1000, beta = 1, sigma=1, keep.all=FALSE, B=B, model="halfnorm"){
    # Function simulates coordinates of individuals on a square
    # Square is [0,2B] x [0,2B], with a count location on the point (B, B)
    #   N: total population size in the square
    #   beta: coefficient of SOEMTHING on spatial covariate x
    #   sigma: scale of half-normal detection function
    #   B: circle radius
    #   keep.all: return the data for y=0 individuals or not
    library(raster)      # Load required packages
    library(plotrix)
    
    # Create coordinates for 30 x 30 grid
    delta <- (2*B-0)/30                # '2D bin width'
    grx <- seq(delta/2, 2*B - delta/2, delta) # mid-point coordinates
    gr <- expand.grid(grx,grx)         # Create grid coordinates
    
    # Create spatially correlated covariate x and plot it
    V <- exp(-e2dist(gr,gr)/1)
    x <- t(chol(V))%*%rnorm(900)
    par(mar=c(3,3,3,6))
    image(rasterFromXYZ(cbind(gr,x)), col=topo.colors(10))
    draw.circle(3, 3, B)
    points(3, 3, pch="+", cex=3)
    
    # Simulate point locations as function of habitat covariate x
    probs <- exp(beta*x)/sum(exp(beta*x)) # probability of point in pixel (sum = 1)
    pixel.id <- sample(1:900, N, replace=TRUE, prob=probs)
    # could simulate randomly within the pixel but it won't matter so place centrally
    u1 <- gr[pixel.id,1]
    u2 <- gr[pixel.id,2]
    points(u1, u2, pch=20, col='black', cex = 0.8)  # plot points
    title("This is so cool !")         # express your appreciation of all this
    
    d <- sqrt((u1 - B)^2 + (u2-B)^2)   # distance to center point of square
    #plot(u1, u2, pch = 1, main = "Point transect")
    N.real <- sum(d<= B)               # Population size inside of count circle
    
    # Can only count individuals in the circle, so set to zero detection probability of individuals in the corners (thereby truncating them)
    # p <- ifelse(d< B, 1, 0) * exp(-d*d/(2*(sigma^2)))
    # We do away with the circle constraint here.   
    if(model=="hazard")
      p <- 1-exp(-exp(-d*d/(2*sigma*sigma)))
    if(model=="halfnorm")
      p <- exp(-d*d/(2*sigma*sigma))
    # Now we decide whether each individual is detected or not
    y <- rbinom(N, 1, p)                                           # detected or not
    points(u1[d<= B], u2[d<= B], pch = 16, col = "black", cex = 1) # not detected
    points(u1[y==1], u2[y==1], pch = 16, col = "red", cex = 1)     # detected
    
    # Put all of the data in a matrix
    if(!keep.all){
      u1 <- u1[y==1]
      u2 <- u2[y==1]
      d <- d[y==1]   }
    # Output
    return(list(model=model, N=N, beta=beta, B=B, u1=u1, u2=u2, d=d, y=y, N.real=N.real, Habitat=x, grid=gr))
  }

# Generate one data set and harvest the output
set.seed(1234)
str(tmp <- sim.spatialDS(N=200, beta=1, sigma=1.5, keep.all=FALSE, B=3)) # Fig. 9-7


# Harvest data
B <- tmp$B
d <- tmp$d
u1 <- tmp$u1
u2 <- tmp$u2
Habitat <- as.vector(tmp$Habitat)
Habitat <- Habitat - mean(Habitat)
Habgrid <- tmp$grid
nind <- length(d)
G <- nrow(Habgrid)

# Do data augmentation, including for pixel ID
M <- 400
nz <- M-nind
pixel <- rep(NA, M)   # We use discrete "pixel ID" here instead of "s"
y <- c(rep(1,nind), rep(0,nz))

# Pick some starting values and figure out the pixel of each observation
s <- cbind(u1,u2)
s <- rbind(s, matrix(NA,nrow=nz,ncol=2))
D <- e2dist(s[1:nind,], Habgrid)
for(i in 1:nind){
  pixel[i] <- (1:ncol(D))[D[i,]==min(D[i,])]
}

# Bundle and summarize the data for BUGS
str(data <- list (B=B, 
                  nind=nind, 
                  y=y, 
                  nz=nz, Habitat=Habitat, Habgrid=Habgrid, G=G, pixel = pixel))

library(nimble)
# Write BUGS model
# "spatialDS.txt" in the book code
Section9p8p3_code <- nimbleCode({ 
    
    # Prior distributions
    sigma ~ dunif(0,10)
    psi ~ dunif(0,1)
    beta ~ dnorm(0,0.01)
    
    for(g in 1:G){   # g is the pixel index, there are G total pixels
      probs.num[g] <- exp(beta*Habitat[g])
      probs[g] <- probs.num[g]/sum(probs.num[1:G])
    }
    
    # Models for DA variables and location (pixel)
    for(i in 1:(nind+nz)){
      z[i] ~ dbern(psi)
      pixel[i] ~ dcat(probs[1:G])
      s[i,1:2] <- Habgrid[pixel[i], ]   # location = derived quantity  
      # compute distance = derived quantity
      d[i] <- pow(pow( s[i,1]-B,2) + pow(s[i,2]-B,2), 0.5)
      p[i] <- exp(-d[i]*d[i]/(2*sigma*sigma))  # Half-normal detetion function
      mu[i]<- p[i]*z[i]
      y[i] ~ dbern(mu[i])                      # Observation model
    }
    # Derived parameters
    N <- sum(z[1:(nind+nz)])           # N is a derived parameter
    D <- N/9                           # area = 9 ha 
    }
)


#  MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nthin <- 2   ;   nc <- 3

# Create inits and define parameters to monitor
inits <- function(){  list (sigma=runif(1,1,10),
                            psi=runif(1),
                            beta=0,
                            z = c(rep(1,nind), 
                                  rep(0,nz)),
                            # Add inits for pixel to avoid nimble warnings
                            pixel = c(rep(NA, nind),
                                      sample(1:G, size = nz, replace = TRUE))) }

params <- c("sigma", "N", "psi", "beta", "D")

# The stochastic indexing of Habgrid requires that we separate constants
# from data
constants <- data[c('B', 'nind', 'nz', 'G')]
data <- data[c('y', 'Habitat', 'Habgrid', 'pixel')]
data$Habgrid <- as.matrix(data$Habgrid)

# Run nimble, check convergence and summarize posteriors
out2 <- nimbleMCMC(
  code = Section9p8p3_code,
  constants = constants,
  data = data,
  inits = inits,
  monitors = params,
  nburnin = 2*nb,
  niter = 2*ni,
  samplesAsCodaMCMC = TRUE
)

## look at posterior information
library(mcmcplots)
colnames(out2)
mcmcplot(out2)

Section9p8p3_compare <- compareMCMCs(
  modelInfo = list(
    code = Section9p8p3_code,
    constants = constants,
    data = data,
    inits = inits()
  ),
  MCMCs = c('nimble', 'jags'),
  monitors = params,
  summary = FALSE,
  burnin = 2*nb,
  niter = 2*ni
)

make_MCMC_comparison_pages(Section9p8p3_compare, modelNames = "Section9p8p3", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section9p8p3.html"))

