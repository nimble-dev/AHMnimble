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
source("Chapter_5_setup.R")
# 5.6 Linear model with normal response (normal GLM): analysis of covariance (ANCOVA)
# -----------------------------------------------------------------------------------

par(mfrow = c(1, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(Cmean ~ factor(facFor), col = c("red", "blue", "green", "grey"), xlab = "Forest cover class", ylab = "Mean count of great tits", frame.plot = F, ylim = c(0,20))
text(0.8, 20, "A", cex=1.6)


## Bundle data
win.data <- list(Cmean = Cmean, M = length(Cmean), elev = elev, facFor = facFor)

require(nimble)
## Specify model in BUGS language for effects parameterizsation:
section5p6_code <- nimbleCode( {
    # Priors
    alpha ~ dnorm(0, 1.0E-06)            # Prior for intercept = effect of level 1 of forest factor
    beta2 ~ dnorm(0, 1.0E-06)            # Prior for slope = effect of elevation for level 1 of forest factor
    beta1[1] <- 0                        # Set to zero effect of first level of facFor
    beta3[1] <- 0                        # Set to zero effect of first level of facFor of elevation
    for(k in 2:4){
    beta1[k] ~ dnorm(0, 1.0E-06)       # Prior for effects of factor facFor
    beta3[k] ~ dnorm(0, 1.0E-06)       # Prior for effects of factor facFor
    }
    tau <- pow(sd, -2)
    sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale
    
    # Likelihood
    for (i in 1:M){
    Cmean[i] ~ dnorm(mu[i], tau)          # precision tau = 1 / variance
    mu[i] <- alpha + beta1[facFor[i]] + beta2 * elev[i] + beta3[facFor[i]] * elev[i]
    }
})

## Initial values
inits <- function() list(alpha = rnorm(1,,10), beta1 = c(NA, rnorm(3,,10)), beta2 = rnorm(1,,10), beta3 = c(NA, rnorm(3,,10)))

## Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "sd")

## MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

## Build, configure and compile model and MCMC:
out3 <- nimbleMCMC(code = section5p6_code,
                     constants = win.data,
                     inits = inits,
                     monitors = params,
                     nburnin = nb,
                     niter = ni)

require(coda)
out3.mcmc <- as.mcmc(out3)
print(out3.mcmc, 3)
require(mcmcplots)
mcmcplot(out3.mcmc)

# Fit model using least-squares (produces MLEs)
(fm <- summary(lm(Cmean ~ as.factor(facFor)*elev)))

# Specify model in BUGS language for means parameterization:
# This is for file "ANCOVA2.txt" in original book code.
section5p6_code.2 <- nimbleCode({
    
    # Priors
    for(k in 1:4){
    alpha[k] ~ dnorm(0, 1.0E-06)       # Priors for intercepts
    beta[k] ~ dnorm(0, 1.0E-06)        # Priors for slopes
    }
    tau <- pow(sd, -2)
    sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale
    
    # Likelihood
    for (i in 1:M){
    Cmean[i] ~ dnorm(mu[i], tau)          # precision tau = 1 / variance
    mu[i] <- alpha[facFor[i]] + beta[facFor[i]] * elev[i]
    }
    
    # Derived quantities: comparison of slopes (now you can forget the delta rule !)
    for(k in 1:4){
    diff.vs1[k] <- beta[k] - beta[1]    # Differences relative to beta[1]
    diff.vs2[k] <- beta[k] - beta[2]    # ... relative to beta[2]
    diff.vs3[k] <- beta[k] - beta[3]    # ... relative to beta[3]
    diff.vs4[k] <- beta[k] - beta[4]    # ... relative to beta[4]
  }
})


# Initial values
inits <- function() list(alpha = rnorm(4,,10), beta = rnorm(4,,10), sd = rnorm(1, 5))

# Parameters monitored
params <- c("alpha", "beta", "sd", "diff.vs1", "diff.vs2", "diff.vs3", "diff.vs4")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

## Build, configure and compile model and MCMC:
out4 <- nimbleMCMC(code = section5p6_code.2,
                   constants = win.data,
                   inits = inits,
                   monitors = params,
                   nburnin = nb,
                   niter = ni)

out4.mcmc <- as.mcmc(out4)
mcmcplot(out4.mcmc)

# Fit model using least-squares (produces MLEs)
(fm <- summary(lm(Cmean ~ as.factor(facFor)*elev-1-elev)))

plot(elev[facFor==1], Cmean[facFor==1], col = "red", ylim = c(0, 20), xlab = "Elevation", ylab = "", frame.plot = F)
points(elev[facFor==2], Cmean[facFor==2], col = "blue")
points(elev[facFor==3], Cmean[facFor==3], col = "green")
points(elev[facFor==4], Cmean[facFor==4], col = "black")
abline(fm$coef[1,1], fm$coef[5,1], col = "red")
abline(fm$coef[2,1], fm$coef[6,1], col = "blue")
abline(fm$coef[3,1], fm$coef[7,1], col = "green")
abline(fm$coef[4,1], fm$coef[8,1], col = "black")
text(-0.8, 20, "B", cex=1.6)

## NIMBLE does not make a list in the same format, so the following code
## has been modified from AHM book code.
##attach.bugs(out4)     # Allows to directly address the sims.list
str(out4[,"diff.vs3[1]"])
par(mfrow = c(1, 3), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
hist(out4[,"diff.vs3[1]"], col = "grey", breaks = 100, main = "", freq=F, ylim = c(0, 0.8))
abline(v = 1, lwd = 3, col = "red")
text(-1.2, 0.8, "A", cex = 2)
hist(out4[,"diff.vs3[2]"], col = "grey", breaks = 100, main = "", freq=F, ylim = c(0, 0.8))
abline(v = 1, lwd = 3, col = "red")
text(-1.4, 0.8, "B", cex = 2)
hist(out4[,"diff.vs3[4]"], col = "grey", breaks = 100, main = "", freq=F, ylim = c(0, 0.8))
abline(v = 1, lwd = 3, col = "red")
text(-2.2, 0.8, "C", cex = 2)

# Prob. difference greater than 1
mean(out4[,"diff.vs3[1]"] > 1)
mean(out4[,"diff.vs3[2]"] > 1)
mean(out4[,"diff.vs3[4]"] > 1)

# compare MCMC efficiecy between NIMBlE and jags:

Section5p6_compare <- compareMCMCs(
  modelInfo <- list(
    code = section5p6_code.2,
    data = win.data,
    inits = inits()),
  MCMCs <- c('nimble', 'jags'),
  summary = TRUE,
  burnin = nb,
  niter = ni
)

make_MCMC_comparison_pages(Section5p6_compare, modelNames = "Section5p6", dir = outputDirectory)

browseURL(Sys.glob(file.path(outputDirectory, "Section5p6.html")))
## Again, we know that for simple linear models in particular, JAGS does well,
## but these are not problematic models for MCMC in general.
