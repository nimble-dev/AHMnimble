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
# This file requires that Section_5p3_example_nimble.R has been run.

# 5.5.1 Some responses missing
# ------------------------------------------------------------------------
## Copy mean counts and turn first 10 into NAs
Cm <- Cmean         ## Copy Cmean into Cm
Cm[1:10] <- NA      ## turn first 10 into missing

## Bundle data (inside BUGS use Cm for Cmean)
win.data <- list(Cmean = Cm, M = length(Cm), elev = elev, forest = forest)

params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "Cmean", "mu")

## The following from AHM book code will not work in NIMBLE.  NIMBLE
## does not monitor subsets of a variable.
## params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "Cmean[1:10]", "mu[1:10]")

inits <- function() list(alpha0 = rnorm(1, 0, 10), alpha1 = rnorm(1, 0, 10),    
                         alpha2 = rnorm(1, 0, 10), alpha3 = rnorm(1, 0, 10))

require(nimble)
## Here we will use the nimbleMCMC function rather than write out all of the steps.
out1.1 <- nimbleMCMC(code = section5p3_code,
                     constants = win.data,
                     inits = inits,
                     monitors = params,
                     nburnin = nb,
                     niter = ni)

require(coda)
out1.1.mcmc <- as.mcmc(out1.1) ## comparison code modified for nimble/coda objects
summary(out1.1.mcmc[,1:10])[[2]]
print(cbind(Truth = Cmean[1:10], summary(out1.1.mcmc[,1:10])[[2]]))

