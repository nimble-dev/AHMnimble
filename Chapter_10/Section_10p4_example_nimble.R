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
source("Chapter_10_setup.R")
# 10.4 A slightly more complex site-occupancy model with covariates
# ------------------------------------------------------------------------

# Choose sample sizes and prepare obs. data array y
set.seed(1)                   # So we all get same data set
M <- 100                      # Number of sites
J <- 3                        # Number of presence/absence measurements
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Create a covariate called vegHt
vegHt <- sort(runif(M, -1, 1)) # sort for graphical convenience

# Choose parameter values for occupancy model and compute occupancy
beta0 <- 0                    # Logit-scale intercept
beta1 <- 3                    # Logit-scale slope for vegHt
psi <- plogis(beta0 + beta1 * vegHt) # Occupancy probability
# plot(vegHt, psi, ylim = c(0,1), type = "l", lwd = 3) # Plot psi relationship

# Now visit each site and observe presence/absence perfectly
z <- rbinom(M, 1, psi)        # True presence/absence

# Look at data so far
table(z)

# Plot the true system state
par(mfrow = c(1, 3), mar = c(5,5,2,2), cex.axis = 1.5, cex.lab = 1.5)
plot(vegHt, z, xlab="Vegetation height", ylab="True presence/absence (z)", frame = F, cex = 1.5)
plot(function(x) plogis(beta0 + beta1*x), -1, 1, add=T, lwd=3, col = "red")

# Create a covariate called wind
wind <- array(runif(M * J, -1, 1), dim = c(M, J))

# Choose parameter values for measurement error model and compute detectability
alpha0 <- -2                        # Logit-scale intercept
alpha1 <- -3                        # Logit-scale slope for wind
p <- plogis(alpha0 + alpha1 * wind) # Detection probability
# plot(p ~ wind, ylim = c(0,1))     # Look at relationship

# Take J = 3 presence/absence measurements at each site
for(j in 1:J) {
  y[,j] <- rbinom(M, z, p[,j])
}
sum(apply(y, 1, max))               # Number of sites with observed presences

# Plot observed data and true effect of wind on detection probability
plot(wind, y, xlab="Wind", ylab="Observed det./nondetection data (y)", frame = F, cex = 1.5)
plot(function(x) plogis(alpha0 + alpha1*x), -1, 1, add=T, lwd=3, col = "red")

# Look at the data: occupancy, true presence/absence (z), and measurements (y)
cbind(psi=round(psi,2), z=z, y1=y[,1], y2=y[,2], y3=y[,3])

# Create factors
time <- matrix(rep(as.character(1:J), M), ncol = J, byrow = TRUE)
hab <- c(rep("A", 33), rep("B", 33), rep("C", 34))  # Must have M = 100

# Bundle and summarize data set
str( win.data <- list(y = y, vegHt = vegHt, wind = wind, M = nrow(y), J = ncol(y), XvegHt = seq(-1, 1, length.out=100), Xwind = seq(-1, 1, length.out=100)) )

library(nimble)
# Specify model in BUGS language
# "model.txt" in book code
Section10p4_code <- nimbleCode({
    
    # Priors
    mean.p ~ dunif(0, 1)         # Detection intercept on prob. scale
    alpha0 <- logit(mean.p)      # Detection intercept
    alpha1 ~ dunif(-20, 20)      # Detection slope on wind
    mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
    beta0 <- logit(mean.psi)     # Occupancy intercept
    beta1 ~ dunif(-20, 20)       # Occupancy slope on vegHt
    
    # Likelihood
    for (i in 1:M) {
      # True state model for the partially observed true state
      z[i] ~ dbern(psi[i])      # True occupancy z at site i
      logit(psi[i]) <- beta0 + beta1 * vegHt[i]
      for (j in 1:J) {
        # Observation model for the actual observations
        y[i,j] ~ dbern(p.eff[i,j])    # Detection-nondetection at i and j
        p.eff[i,j] <- z[i] * p[i,j]   # 'straw man' for WinBUGS
        logit(p[i,j]) <- alpha0 + alpha1 * wind[i,j]
      }
    }
    
    # Derived quantities
    N.occ <- sum(z[1:M])       # Number of occupied sites among sample of M
    psi.fs <- N.occ/M       # Proportion of occupied sites among sample of M
    for(k in 1:100){
      logit(psi.pred[k]) <- beta0 + beta1 * XvegHt[k] # psi predictions
      logit(p.pred[k]) <- alpha0 + alpha1 * Xwind[k]  # p predictions
    }
  }
)


# Initial values: must give for same quantities as priors given !
zst <- apply(y, 1, max)        # Avoid data/model/inits conflict
inits <- function(){list(z = zst, 
                         mean.p = runif(1), 
                         alpha1 = runif(1), 
                         mean.psi = runif(1), 
                         beta1 = runif(1))}

# Parameters monitored
params <- c("alpha0", "alpha1", "beta0", "beta1", "N.occ", "psi.fs", "psi.pred", "p.pred", "z") # Also estimate z = "conditional occ. prob."

# MCMC settings
ni <- 25000   ;   nt <- 10   ;   nb <- 2000   ;   nc <- 3

# Call nimble from R and summarize posteriors
out1B <- nimbleMCMC(
  code = Section10p4_code,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)

library(mcmcplots)
colnames(out1B)
params_of_interest <- c("alpha0", "alpha1", "beta0", "beta1", "N.occ", "psi.fs")
mcmcplot(out1B[, params_of_interest])
## Mixing looks good

Section10p4_compare <- compareMCMCs(
  modelInfo = list(
    code = Section10p4_code,
    data = win.data,
    inits = inits()
  ),
  ## making the comparison for all latent states in params would be cumbersome
  monitors = params_of_interest,
  MCMCs = c('nimble', 'jags'),
  summary = FALSE,
  burnin = nb,
  niter = ni
)

make_MCMC_comparison_pages(Section10p4_compare, modelNames = "Section10p4", dir = outputDirectory)

browseURL(Sys.glob(file.path(outputDirectory, "Section10p4.html")))
