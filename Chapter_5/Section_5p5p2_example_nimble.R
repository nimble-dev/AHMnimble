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

## All responses are missing:

## Bundle data: simply drop the response from list
win.data <- list(M = length(Cm), elev = elev, forest = forest)

## Alternatively, add all-NA data vector
win.data$Cmean <- as.numeric(rep(NA, length(Cmean)))
str(win.data)  # Cmean is numeric

## Parameters monitored
##params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "mu[1:2]")
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "mu")

## Take steps to build, configure and compile MCMC.

out1.2 <- nimbleMCMC(code = section5p3_code,
                     constants = win.data,
                     inits = inits,
                     monitors = params,
                     nburnin = nb,
                     niter = ni)

library(coda)
out1.2.mcmc <- as.mcmc(out1.2)

library(mcmcplots)
mcmcplot(out1.2.mcmc[,1:6])
## These chains are not very well mixed.
## They could be run longer, or different samplers could be used.

print(out1.2.mcmc[,1:6], 2)

## Further plots from AHM book code do not work due to different MCMC object formats.

