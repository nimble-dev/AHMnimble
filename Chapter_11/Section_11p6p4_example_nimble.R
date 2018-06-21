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
# NIMBLE has not been checked for this case.
# 11.6.4 Modeling species richness in a two-step analysis
# ------------------------------------------------------------------------
# Extract estimates of N from model 5
N.pm <- out5$summary[291:557, 1]       # Posterior means of Nsite
N.pm <- out5$Nsite
N.psd <- out5$summary[291:557, 2]      # ... posterior sd's of Nsite
N.cri <- out5$summary[291:557, c(3,7)] # ... CRL's of Nsite

# Plot estimates as a function of elevation
elev <- data$elev[1:267]
plot(elev, N.pm, xlab = "Altitude (m a.s.l.)", ylab = "Estimated avian species richness", ylim = c(0, 70), frame = F)
segments(elev, N.cri[,1], elev, N.cri[,2], col = "grey")
lines(smooth.spline(N.pm ~ elev, w = 1 / N.psd), col = "grey", lwd = 3)


# Bundle and summarize data set
pred.ele <- (seq(200, 2750,5) - mean.ele) / sd.ele # elevation standardised
str(win.data <- list(ele = ele, N = N.pm, psd = N.psd, n = length(N.pm), pred.ele = pred.ele, npred = length(pred.ele)))

# Define model in BUGS language
Section11p6p4_code <- nimbleCode({
    
    # Priors 
    for(v in 1:4){         # Priors for intercept and polynomial coefficients
      beta[v] ~ dnorm(0, 0.0001)
    }
    tau.site <- pow(sd.site, -2)
    sd.site ~ dunif(0,10)
    
    # Likelihood
    for (i in 1:n){ 
      N[i] ~ dnorm(muN[i], tau.psd[i]) # Measurement error model for estimated N
      tau.psd[i] <- pow(psd[i], -2)    # 'Known' part of residual: meas. error
      muN[i] <- beta[1] + beta[2] * ele[i] + beta[3] * pow(ele[i],2) + 
      beta[4] * pow(ele[i],3) + eps.site[i] # add another source of uncertainty
      eps.site[i] ~ dnorm(0, tau.site) # this is the usual 'residual'
    }
    # Get predictions for plot
    for(i in 1:npred){
      Npred[i] <- beta[1] + beta[2] * pred.ele[i] + beta[3] * pow(pred.ele[i],2) + beta[4] * pow(pred.ele[i],3)
    }
  } # end model 
)

# Initial values, params monitored, and MCMC settings
inits <- function() list(beta = rnorm(4)) 
params <- c("beta", "sd.site", "Npred") 
ni <- 12000   ;   nt <- 10   ;   nb <- 2000   ;   nc <- 3 

# Call nimble and summarize posterior
out <- nimbleMCMC(
  code = Section11p6p4_code,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = 4*ni,
  samplesAsCodaMCMC = TRUE
)

library(coda)
library(mcmcplots)

