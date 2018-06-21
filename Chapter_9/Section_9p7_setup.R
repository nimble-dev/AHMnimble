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
# 9.7.1 Simulating the ISSJ data over multiple years
# ------------------------------------------------------------------------
# We load the ISSJ data analyzed in chapter 8, package into an unmarked frame  
library(unmarked)
library(AHMbook)
data(issj)
covs <- issj[,c("elevation","forest","chaparral")]
area <- pi*300^2 / 100^2             # Area in ha
jayumf <- unmarkedFrameGDS(y=as.matrix(issj[,1:3]),
                           siteCovs=data.frame(covs, area), numPrimary=1,
                           dist.breaks=c(0, 100, 200, 300),
                           unitsIn="m", survey="point")
sc <- siteCovs(jayumf)
sc.s <- scale(sc)
sc.s[,"area"] <- pi*300^2 / 10000  # Don't standardize area
covs<- siteCovs(jayumf) <- sc.s
summary(jayumf)

# Fit the model using gdistsamp and look at the fit summary
(nb.C2E.C <- gdistsamp( ~chaparral + I(chaparral^2) + elevation , ~1, ~chaparral,
                        data =jayumf, output="abund", mixture="NB", K = 150))

# Get coefficient estimates to be used in data simulation
beta <- coef(nb.C2E.C)
betaFall <- beta[c("lambda(Int)", "lambda(chaparral)",
                   "lambda(elevation)", "lambda(I(chaparral^2))")]

# Predict expected abundance per point count on log-scale for simulation
Xmat <- cbind(rep(1,307),covs[,3],covs[,3]^2,covs[,1]) # Order: chap, chap^2, elev
loglam <- Xmat%*%(betaFall)
lamnew <- exp(loglam)

# Parameters of the detection function
dparm <- beta[c("p(Int)", "p(chaparral)")]
sigma <- exp(Xmat[,c(1, 2)]%*%dparm)
J <- nsites <- 307 # number of sampling points

# Number of years
nyrs <- 6

# Set dynamics parameters to achieve a target growth rate of 0.95
phi <- 0.6       # Survival probability
gamma <- 0.35    # Recruitment rate

# Distance category info
db <- c(0,50, 100, 150, 200, 250, 300)
midpt <- c(25, 75, 125, 175, 225, 275)
nD <- length(midpt)
delta <- 50       # Distance interval width
B <- 300

# Simulate an ISSJ data set and harvest the data objects
set.seed(2015)
dat <- issj.sim(B=300, db = db, lam=lamnew, sigma=sigma, phi=phi, gamma=gamma, npoints=nsites, nyrs=nyrs)

y <- dat$y
dclass <- dat$dclass
site <- dat$site

# Bundle and summarize the data set
str(data1<-list(nsites=nsites, chap=as.vector(covs[,"chaparral"])[dat$cell], 
                chap2=as.vector(covs[,"chaparral"]^2)[dat$cell],
                elev=as.vector(covs[,"elevation"])[dat$cell], T=nyrs, nD=nD, midpt=midpt, 
                B=B, delta=delta, y=y, dclass=dclass, site=site, nind=sum(y)) )



