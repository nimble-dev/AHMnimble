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
# 9.2.1 Simulating HDS data with group size
# ------------------------------------------------------------------------
# This file contains code which creates data which is necessary other examples in this chapter

# Function to simulate data under HDS protocol with groups
simHDSg <- function(type = "line", nsites = 100, lambda.group = 0.75, alpha0 = 0, alpha1 = 0.5, beta0 = 1, beta1 = 0.5, B = 4, discard0 = TRUE){
  #
  # Function simulates hierarchical distance sampling (HDS) data for groups under 
  #   either a line (type = "line") or a point (type = "point") transect protocol 
  #   and using a half-normal detection function (Buckland et al. 2001).
  #   Other function arguments:
  #     nsites: Number of sites (spatial replication)
  #     lambda.group: Poisson mean of group size
  #     alpha0, alpha1: intercept and slope of log-linear model relating sigma of
  #        half-normal detection function to group size
  #     beta0, beta1: intercept and slope of log-linear model relating the Poisson
  #        mean of the number of groups per unit area to habitat
  #     B: strip half width
  #
  # Get covariates
  habitat <- rnorm(nsites)              # Simulated covariate  
  
  # Simulate abundance model for groups (Poisson GLM for N)
  lambda <- exp(beta0 + beta1*habitat)  # Density of groups per "square"
  N <- rpois(nsites, lambda)            # site-specific number of groups
  N.true <- N                           # for point: inside of B
  
  # Simulate observation model
  data <- groupsize <- NULL
  
  for(i in 1:nsites){
    if(N[i]==0){
      data <- rbind(data,c(i,NA,NA,NA,NA,NA)) # save site, y=1, u, v, d
      next
    }
    
    if(type=="line"){
      # Simulation of distances, uniformly, for each individual in the population
      d <- runif(N[i], 0, B)
      gs <- rpois(N[i],lambda.group) +1  # Observable group sizes >= 1
      groupsize<-c(groupsize,gs)
      sigma.vec <- exp(alpha0 + alpha1*(gs-1))  # Subtract 1 for interpretation
      # Detection probability for each group
      p <- exp(-d*d/(2*(sigma.vec^2)))
      # Determine if individuals are captured or not
      y <- rbinom(N[i], 1, p)
      u1 <- u2 <- rep(NA,N[i])
      # Subset to "captured" individuals only
      d <- d[y==1] ; u1 <- u1[y==1] ; u2 <- u2[y==1] ; gs <- gs[y==1] ; y <- y[y==1]
    }
    
    if(type=="point"){
      # Simulation of data on a circle of radius B (algorithm of Wallin)
      angle <- runif(N[i], 0, 360)
      r2 <- runif(N[i], 0, 1)
      r <-  B*sqrt(r2)
      u1 <-  r*cos(angle)  + B
      u2 <-  r*sin(angle)  + B
      
      d <- sqrt((u1 - B)^2 + (u2-B)^2)
      N.true[i] <- sum(d<= B)    # Population size inside of count circle, should be N[i] here.
      gs <- rpois(N[i], lambda.group) + 1
      groupsize <-c(groupsize,gs)
      sigma.vec <- exp(alpha0 + alpha1*(gs-1))
      # For counting individuals on a circle so we truncate p here
      p <- ifelse(d<(B), 1, 0)*exp(-d*d/(2*(sigma.vec^2)))
      y <- rbinom(N[i], 1, p)
      # Subset to "captured" individuals only
      d <- d[y==1] ; u1 <- u1[y==1] ; u2 <- u2[y==1] ; gs <- gs[y==1] ; y <- y[y==1]
    }
    # Now compile things into a matrix and insert NA if no individuals were 
    # captured at site i. Coordinates (u,v) are preserved.
    if(sum(y) > 0)
      data <- rbind(data,cbind(rep(i, sum(y)), y, u1, u2, d, gs))
    else
      data <- rbind(data,c(i,NA,NA,NA,NA,NA)) # make a row of missing data
  }
  # Subset to sites at which individuals were captured. You may or may not
  #  do this depending on how the model is formulated so be careful. 
  if(discard0)
    data <- data[!is.na(data[,2]),]
  
  # Visualisation
  if(type=="line"){       # For line transect
    par(mfrow = c(1, 3))
    hist(data[,"d"], col = "lightblue", breaks = 20, main = 
           "Frequency of distances to groups", xlab = "Distance")
    ttt <- table(data[,1])
    n <- rep(0, nsites)
    n[as.numeric(rownames(ttt))] <- ttt
    plot(habitat, n, main = "Observed group counts (n) vs. habitat", frame = F)
    plot(table(data[,"gs"]), main = "Observed group sizes", ylab = "Frequency", frame = F)
  }
  
  if(type=="point"){       # For point transect
    par(mfrow = c(2,2))
    plot(data[,"u1"], data[,"u2"], pch = 16, main = 
           "Located groups in point transects", xlim = c(0, 2*B),
         ylim = c(0, 2*B), col = data[,1], asp = 1)
    points(B, B, pch = "+", cex = 3)
    library(plotrix)
    draw.circle(B, B, B)
    hist(data[,"d"], col = "lightblue", breaks = 20, main = 
           "Frequency of distances to groups", xlab = "Distance")
    ttt <- table(data[,1])
    n <- rep(0, nsites)
    n[as.numeric(rownames(ttt))] <- ttt
    plot(habitat, n, main = "Observed group counts (n) vs. habitat", frame = F)
    plot(table(data[,"gs"]), main = "Observed group sizes", ylab = "Frequency", frame = F)
  }
  
  # Output
  list(type = type, nsites = nsites, lambda.group = lambda.group, alpha0 = alpha0, alpha1 = alpha1, beta0 = beta0, beta1 = beta1, B = B, data=data, habitat=habitat, N = N, N.true = N.true, groupsize=groupsize)
}


data <- simHDSg(type = "line")     # Defaults for line transect data
data <- simHDSg(type = "point")    # Default for point transect data
data <- simHDSg(lambda.group = 5)  # Much larger groups
data <- simHDSg(lambda.group = 5, alpha1 = 0) # No effect of groups size on p
