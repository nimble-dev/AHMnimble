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

source("Chapter_5_setup_nimble.R")

# This file contains code needed to set up examples from Chapter 5.

# 4.3 Packaging everything in a function
# ------------------------------------------------------------------------
# Function definition with set of default values
data.fn <- function(M = 267, J = 3, mean.lambda = 2, beta1 = -2, beta2 = 2, beta3 = 1, mean.detection = 0.3, alpha1 = 1, alpha2 = -3, alpha3 = 0, show.plot = TRUE){
  #
  # Function to simulate point counts replicated at M sites during J occasions.
  # Population closure is assumed for each site.
  # Expected abundance may be affected by elevation (elev),
  # forest cover (forest) and their interaction.
  # Expected detection probability may be affected by elevation,
  # wind speed (wind) and their interaction.
  # Function arguments:
  #     M: Number of spatial replicates (sites)
  #     J: Number of temporal replicates (occasions)
  #     mean.lambda: Mean abundance at value 0 of abundance covariates
  #     beta1: Main effect of elevation on abundance
  #     beta2: Main effect of forest cover on abundance
  #     beta3: Interaction effect on abundance of elevation and forest cover
  #     mean.detection: Mean detection prob. at value 0 of detection covariates
  #     alpha1: Main effect of elevation on detection probability
  #     alpha2: Main effect of wind speed on detection probability
  #     alpha3: Interaction effect on detection of elevation and wind speed
  #     show.plot: if TRUE, plots of the data will be displayed;
  #        set to FALSE if you are running simulations.
  
  
  # Create covariates
  elev <- runif(n = M, -1, 1)                         # Scaled elevation
  forest <- runif(n = M, -1, 1)                       # Scaled forest cover
  wind <- array(runif(n = M*J, -1, 1), dim = c(M, J)) # Scaled wind speed
  
  # Model for abundance
  beta0 <- log(mean.lambda)               # Mean abundance on link scale
  lambda <- exp(beta0 + beta1*elev + beta2*forest + beta3*elev*forest)
  N <- rpois(n = M, lambda = lambda)      # Realised abundance
  Ntotal <- sum(N)                        # Total abundance (all sites)
  psi.true <- mean(N>0)                   # True occupancy in sample
  
  # Plots
  if(show.plot){
    par(mfrow = c(2, 2), cex.main = 1)
    devAskNewPage(ask = TRUE)
    curve(exp(beta0 + beta1*x), -1, 1, col = "red", main = "Relationship lambda-elevation \nat average forest cover", frame.plot = F, xlab = "Scaled elevation")
    plot(elev, lambda, xlab = "Scaled elevation", main = "Relationship lambda-elevation \nat observed forest cover", frame.plot = F)
    curve(exp(beta0 + beta2*x), -1, 1, col = "red", main = "Relationship lambda-forest \ncover at average elevation", xlab = "Scaled forest cover", frame.plot = F)
    plot(forest, lambda, xlab = "Scaled forest cover", main = "Relationship lambda-forest cover \nat observed elevation", frame.plot = F)
  }
  
  # Model for observations
  alpha0 <- qlogis(mean.detection)        # mean detection on link scale
  p <- plogis(alpha0 + alpha1*elev + alpha2*wind + alpha3*elev*wind)
  C <- matrix(NA, nrow = M, ncol = J)     # Prepare matrix for counts
  for (i in 1:J){                         # Generate counts by survey
    C[,i] <- rbinom(n = M, size = N, prob = p[,i])
  }
  summaxC <- sum(apply(C,1,max))          # Sum of max counts (all sites)
  psi.obs <- mean(apply(C,1,max)>0)       # Observed occupancy in sample
  
  # More plots
  if(show.plot){
    par(mfrow = c(2, 2))
    curve(plogis(alpha0 + alpha1*x), -1, 1, col = "red", main = "Relationship p-elevation \nat average wind speed", xlab = "Scaled elevation", frame.plot = F)
    matplot(elev, p, xlab = "Scaled elevation", main = "Relationship p-elevation\n at observed wind speed", pch = "*", frame.plot = F)
    curve(plogis(alpha0 + alpha2*x), -1, 1, col = "red", main = "Relationship p-wind speed \n at average elevation", xlab = "Scaled wind speed", frame.plot = F)
    matplot(wind, p, xlab = "Scaled wind speed", main = "Relationship p-wind speed \nat observed elevation", pch = "*", frame.plot = F)
    
    matplot(elev, C, xlab = "Scaled elevation", main = "Relationship counts and elevation", pch = "*", frame.plot = F)
    matplot(forest, C, xlab = "Scaled forest cover", main = "Relationship counts and forest cover", pch = "*", frame.plot = F)
    matplot(wind, C, xlab = "Scaled wind speed", main = "Relationship counts and wind speed", pch = "*", frame.plot = F)
    desc <- paste('Counts at', M, 'sites during', J, 'surveys')
    hist(C, main = desc, breaks = 50, col = "grey")
  }
  
  # Output
  return(list(M = M, J = J, mean.lambda = mean.lambda, beta0 = beta0, beta1 = beta1, beta2 = beta2, beta3 = beta3, mean.detection = mean.detection, alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, elev = elev, forest = forest, wind = wind, lambda = lambda, N = N, p = p, C = C, Ntotal = Ntotal, psi.true = psi.true, summaxC = summaxC, psi.obs = psi.obs))
}

# 5.3 Linear model with normal response (normal GLM): multiple linear regression
# ------------------------------------------------------------------------

# Generate data with data.fn from chapter 4
set.seed(24)
data <- data.fn(show.plot = FALSE) ## added show.plot = FALSE
str(data)
attach(data)

# Summarize data by taking mean at each site and plot
Cmean <- apply(C, 1, mean)
# Following figures from book code are skipped for nimble setup
# par(mfrow = c(1,3))
# hist(Cmean, 50)               # Very skewed
# plot(elev, Cmean)
# plot(forest, Cmean)

# Package the data needed in a bundle
win.data <- list(Cmean = Cmean, M = length(Cmean), elev = elev, forest = forest)
str(win.data)                    # Check what's in win.data

                                       #
# 5.6 Linear model with normal response (normal GLM): analysis of covariance (ANCOVA)
# -----------------------------------------------------------------------------------
# Generate factor and plot raw data in boxplot as function of factor A
facFor <- as.numeric(forest < -0.5)         # Factor level 1
facFor[forest < 0 & forest > -0.5] <- 2     # Factor level 2
facFor[forest < 0.5 & forest > 0] <- 3      # Factor level 3
facFor[forest > 0.5] <- 4                   # Factor level 4
#table(facFor)                               # every site assigned a level OK

# 5.9 Poisson generalized linear model (Poisson GLM)
# --------------------------------------------------

# Summarize data by taking max at each site
Cmax <- apply(C, 1, max)
#table(Cmax)

# 5.11 Binomial generalised linear model (binomial GLM, logistic regression)
# --------------------------------------------------------------------------

# Quantize counts from first survey and describe
y1 <- as.numeric(C[,1] > 0)  # Gets 1 if first count greater than zero

# 5.14 Random-effects binomial GLM (binomial GLMM)
# ------------------------------------------------

# Get detection/nondetection response
y <- C
y[y > 0] <- 1
