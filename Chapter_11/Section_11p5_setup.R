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
source("Section_11p3.R")
# 11.5.1 Simple Poisson regression for the observed community size
# ------------------------------------------------------------------------
# Get covariates and standardise them
# Quadrat elevation and forest cover
data <- MHB2014
orig.ele <- data$sites$elev[1:nsite]
(mean.ele <- mean(orig.ele, na.rm = TRUE))
(sd.ele <- sd(orig.ele, na.rm = TRUE))
ele <- (orig.ele - mean.ele) / sd.ele
orig.forest <- data$sites$forest[1:nsite]
(mean.forest <- mean(orig.forest, na.rm = TRUE))
(sd.forest <- sd(orig.forest, na.rm = TRUE))
forest <- (orig.forest - mean.forest) / sd.forest

# Average date and duration of survey
tmp <- cbind(data$date)[1:nsite,] 
orig.mdate <- apply(tmp, 1, mean, na.rm = TRUE)
(mean.mdate <- mean(orig.mdate[-NAsites]))   # drop unsurved site
(sd.mdate <- sd(orig.mdate[-NAsites]))
mdate <- (orig.mdate - mean.mdate) / sd.mdate
mdate[NAsites] <- 0                 # impute mean for missing

tmp <- cbind(data$dur)[1:nsite,] 
orig.mdur <- apply(tmp, 1, mean, na.rm = TRUE)
(mean.mdur <- mean(orig.mdur[-NAsites]))
(sd.mdur <- sd(orig.mdur[-NAsites]))
mdur <- (orig.mdur - mean.mdur) / sd.mdur
mdur[NAsites] <- 0                  # impute mean for missing

# Get observed species richness per site and rep and plot
CC <- apply(y, c(1,2), sum, na.rm = TRUE)
CC[CC == 0] <- NA            # 0 means not surveyed
matplot(t(CC), type = 'l', lty = 1, lwd = 2, xlab = "First to third survey", ylab = "Number of species detected", frame = F)  # Fig. 11�6 right

# Get survey date and survey duration and standardise both
# Survey date (this is Julian date, with day 1 being April 1)
orig.DAT <- cbind(data$date)[1:nsite,]
(mean.date <- mean(orig.DAT, na.rm = TRUE))
(sd.date <- sd(c(orig.DAT), na.rm = TRUE))
DAT <- (orig.DAT - mean.date) / sd.date      # scale
DAT[is.na(DAT)] <- 0                         # impute missings
# Survey duration (in minutes)
orig.DUR <- cbind(data$dur)[1:nsite,]
(mean.dur <- mean(orig.DUR, na.rm = TRUE))
(sd.dur <- sd(c(orig.DUR), na.rm = TRUE))
DUR <- (orig.DUR - mean.dur) / sd.dur        # scale
DUR[is.na(DUR)] <- 0                         # mean impute missings
