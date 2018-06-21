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
# This file contains some customization of Section 6.4 for NIMBLE.
#
# Before running code in this file, it is necessary to run Section_6p4_example_nimble.R
# at least including creation of win.data, inits, and Section6p4_code.
source("Section_6p4_setup.R")

## This function assigns Automated Factor Slice Samplers to each pair of
## (beta0[k], beta1[k]) and each pair of (alpha0[k], alpha1[k])
makeCustomMCMCconf <- function(model) {
    MCMCconf <- configureMCMC(model)
    MCMCconf$removeSamplers('alpha0')
    MCMCconf$removeSamplers('alpha1')
    MCMCconf$removeSamplers('beta0')
    MCMCconf$removeSamplers('beta1')
    MCMCconf$addSampler(target = c('beta0[1]', 'beta1[1]'), type = "AF_slice")
    MCMCconf$addSampler(target = c('beta0[2]', 'beta1[2]'), type = "AF_slice")
    MCMCconf$addSampler(target = c('beta0[3]', 'beta1[3]'), type = "AF_slice")
    MCMCconf$addSampler(target = c('alpha0[1]', 'alpha1[1]'), type = "AF_slice")
    MCMCconf$addSampler(target = c('alpha0[2]', 'alpha1[2]'), type = "AF_slice")
    MCMCconf$addSampler(target = c('alpha0[3]', 'alpha1[3]'), type = "AF_slice")
    MCMCconf
}

## This compares four MCMCs:
## JAGS
## NIMBLE with default samplers
## NIMBLE with slice samplers
## NIMBLE with the customized samplers described above
 Section6p4_compare <- compareMCMCs(
  modelInfo = list(
    code = Section6p4_code,
    data = win.data,
    inits = inits()
  ),
  MCMCs = c('nimble', 'jags', 'nimble_slice', 'nimble_AFSS'),
  MCMCdefs = list(
      nimble_AFSS = quote({
          makeCustomMCMCconf(Rmodel)
      })),
  summary = TRUE,
  burnin = nb,
  niter = ni,
  monitors = c('alpha0', 'alpha1', 'beta0', 'beta1', 'Nocc', 'Ntotal', 'Nhab', 'N.critical')
)

make_MCMC_comparison_pages(Section6p4_compare, modelNames = "Section6p4", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section6p4.html"))
