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

# 8.2.3 Simulating distance sampling data
# ------------------------------------------------------------------------
strip.width <- 100  # one side of the transect, really half-width  
sigma <- 30         # Scale parameter of half-normal detection function

# Define half-normal detection function
g <- function(x, sig) exp(-x^2/(2*sig^2)) # Function definition
g(30, sig=sigma)    # Detection probability at a distance of 30m

# Plot the detection function
par(mfrow=c(1,2))
curve(g(x, sig=30), 0, 100, xlab="Distance (x)", ylab="Detection prob.", lwd = 2, frame = F)
curve(g(x, sig=60), 0, 100, add=TRUE, lty = 2, lwd = 2)

# Define function to simulate non-hierarchical line transect data
sim.ldata <- function(N = 200, sigma = 30){
  # Function to simulate line transect data under CDS.
  # Function arguments:
  #    N: number of individuals along transect with distance u(-100, 100)
  #    sigma: scale parameter of half-normal detection function
  # Function subjects N individuals to sampling, and then retains the value 
  # of x=distance only for individuals that are captured
  par(mfrow = c(1,2))
  # Plot the detection function
  curve(exp(-x^2/(2*sigma^2)), 0, 100, xlab="Distance (x)", ylab="Detection prob.", lwd = 2, main = "Detection function", ylim = c(0,1))
  text(80, 0.9, paste("sigma:", sigma))
  xall <- runif(N, -100,100) # Distances of all N individuals
  hist(abs(xall), nclass=10, xlab = "Distance (x)", col = "grey", main = "True (grey) \nand observed distances (blue)")
  g <- function(x, sig) exp(-x^2/(2*sig^2))
  p <- g(xall, sig=sigma) # detection probability
  y <- rbinom(N, 1, p) # some inds. are detected and their distance measured
  x <- xall[y==1]      # this has direction (right or left transect side)
  x <- abs(x)          # now it doesn't have direction
  hist(x, col = "blue", add = TRUE)
  return(list(N = N, sigma = sigma, xall = xall, x = x))
}

# Obtain a data set for analysis
set.seed(2015)               # If you want to get same results
tmp <- sim.ldata(sigma = 30) # Execute function and assign results to 'tmp'
attach(tmp)


# Conditional likelihood
Lcond <- function(lsigma){  # Define conditional nll
  sigma<- exp(lsigma)
  -1*sum(log(g(x,sig=sigma)/integrate(g, 0, 100, sig=sigma)$value/100))
}

# Full likelihood
Lfull <- function(parm){    # Define full nll
  sigma <- exp(parm[1])
  n0 <- exp(parm[2])
  N <- length(x)+ n0
  pbar <- integrate(g, 0, 100, sig=sigma)$value/100
  -1*( lgamma(N+1) - lgamma(n0+1) + sum(log(g(x,sig=sigma)/100)) + n0*log(1-pbar) )
}

# Call optim to maximize full likelihood
optim(c(log(30), log(4)), Lfull, hessian=TRUE)


pbar<- integrate(g, 0, 100, sig=exp(3.26))$value/100
n<- length(tmp$x)

(Nhat.condl<- n/pbar)
(Dhat.condl<- Nhat.condl/(10*.2))

n0hat<- exp(5.01)
(Nhat.full<- n + n0hat)
(Dhat.full<- Nhat.full/(10*.2))


# 8.2.4.2 Simulating binned distance sampling data
# ------------------------------------------------------------------------
set.seed(2015)
# Design settings and truth (population size N and detection function g)
interval.width <- 10
strip.width <- 100    # half-width really (one side of transect)
nbins <- strip.width%/%interval.width  
sigma <- 30           # Scale parameter of half-normal detection function
g <- function(x, sig) exp(-x^2/(2*sig^2)) # Half-normal detection function
N <- 200              # Population size

# Method 1: simulate continuous distances and put into intervals
x <- runif(N, -strip.width, strip.width) # Distance all animals
p <- g(x, sig=sigma)  # Detection probability
y <- rbinom(N, 1, p)  # only individuals with y=1 are detected
x <- x[y==1]          # this has direction (right or left side of transect)
x <- abs(x)           # now it doesn't have direction

# Compute the distance category of each observation
xbin <- x %/% interval.width + 1   # note integer division function %/%

# Multinomial frequencies, may have missing levels
y.obs <- table(xbin)

# Pad the frequencies to include those with 0 detections
y.padded <- rep(0,nbins)
names(y.padded) <- 1:nbins
y.padded[names(y.obs)] <- y.obs
y.obs <- y.padded
y.true <- c(y.obs, N-length(xbin)) # Last category is "Not detected"

# Relative frequencies by binning continuous data (pi). These should compare
#  with the cell probabilities computed below when N is very large
(y.rel <- y.true/N)      # Last category is pi(0) from above
(pi0.v1 <- y.rel[nbins+1])

# Compute detection probability in each distance interval
dist.breaks <- seq(0, strip.width, by=interval.width)
p <- rep(NA, length(dist.breaks)-1)
for(j in 1:length(p)){
  p[j] <- integrate(g, dist.breaks[j], dist.breaks[j+1],
                    sig=sigma)$value / (dist.breaks[j+1]-dist.breaks[j]) 
}
round(p, 2)

# Compute the multinomial cell probabilities analytically. These are exact.
# psi = probability of occurring in each interval
interval.width <- diff(dist.breaks)
psi <- interval.width/strip.width
pi <- p * psi
sum(pi)                 # This is 1 � pi(0) from above
(pi0.exact <- 1-sum(pi))

# Method 2: Use rmultinom to simulate binned observations directly
# This includes 0 cells AND n0
pi[length(p)+1] <- 1 - sum(pi)
(y.obs2 <- as.vector(rmultinom(1, N, prob=pi)))
(y.obs2 <- y.obs2[1:nbins]) # Discard last cell for n0 (because not observed)

Lik.binned <- function(parm, data, dist.breaks){
  # Note that the parameters are parm[1] = log(sigma), parm[2] = log(n0)
  
  sigma <- exp(parm[1])
  n0 <- exp(parm[2])
  p <- rep(NA, length(dist.breaks)-1)
  for(j in 1:length(p)) {
    p[j] <- integrate(g, dist.breaks[j], dist.breaks[j+1],
                      sig=sigma)$value / (dist.breaks[j+1]-dist.breaks[j])
  }
  psi <- interval.width/strip.width
  pi <- p * psi
  pi0 <- 1-sum(pi)
  
  N <- sum(data) + n0
  -1*(lgamma(N+1)-lgamma(n0+1) + sum(c(data,n0)*log(c(pi,pi0))))
}

# Evaluate likelihood for some particular value of the parameters
Lik.binned(c(2,0), data=y.obs, dist.breaks=dist.breaks)

# Obtain the MLEs for the simulated data
optim(c(2,0), Lik.binned, data=y.obs, dist.breaks=dist.breaks)



# 8.2.5 Point transect data
# ------------------------------------------------------------------------
# Define function to compute cell probs for binned distance sampling
cp.ri <-function(radius1, radius2, sigma){
  Pi <- 3.141593
  a <- Pi*radius2^2 - Pi*radius1^2
  integrate(function(r, s=sigma) exp(-r^2 / (2 * s^2)) * r, radius1,radius2)$value *(2*Pi/a)  
}

# Define distance intervals and compute multinomial probabilities
delta <- 0.5                   # Width of distance bins
B <- 3                         # Max count distance
dist.breaks <-seq(0, B, delta) # Make the interval cut points
nD <-length(dist.breaks)-1
sigma <- 1
p.x <-rep(NA,nD)               # Conditional detection probabilities
for(i in 1:nD){
  p.x[i] <- cp.ri(dist.breaks[i], dist.breaks[i+1], sigma =1)
}
area <- 3.141593 * dist.breaks[-1]^2
ring.area <- diff(c(0, area))
# Pr(detection| in ring)*Pr(in ring)
cp <- p.x* ring.area/sum(ring.area)



# 8.2.5.1 Simulating point transect data
# ------------------------------------------------------------------------
sim.pdata <- function(N=1000, sigma=1, B=3, keep.all=FALSE) {
  # Function simulates coordinates of individuals on a square
  # Square is [0,2*B] x[0,2*B], with a count location on the center 
  # point (B,B)
  # Function arguments:
  #    N: total population size in the square
  #    sigma: scale of half-normal detection function
  #    B: circle radias
  #    keep.all: return the data for y = 0 individuals or not
  
  # Plot the detection function
  par(mfrow = c(1,2))
  curve(exp(-x^2/(2*sigma^2)), 0, B, xlab="Distance (x)", ylab="Detection prob.", lwd = 2, main = "Detection function", ylim = c(0,1))
  text(0.8*B, 0.9, paste("sigma:", sigma))
  
  # Simulate and plot simulated data
  library(plotrix)
  u1 <-runif(N, 0, 2*B)           # (u1,u2) coordinates of N individuals
  u2 <- runif(N, 0, 2*B)
  d <- sqrt((u1 - B)^2 + (u2 - B)^2) # distance to center point of square
  plot(u1, u2, asp = 1, pch = 1, main = "Point transect")
  N.real <- sum(d<= B)           # Population size inside of count circle
  
  # Can only count indidividuals in the circle, so set to zero detection probability of individuals in the corners (thereby truncating them):
  p <- ifelse(d < B, 1, 0) * exp(-d*d/(2*(sigma^2)))
  # Now we decide whether each individual is detected or not
  y <- rbinom(N, 1, p)
  points(u1[d <= B], u2[d <= B], pch = 16, col = "black")
  points(u1[y==1], u2[y==1], pch = 16, col = "blue")
  points(B, B, pch = "+", cex = 3, col = "red")
  draw.circle(B, B, B)
  
  # Put all of the data in a matrix: 
  #      (note we don't care about y, u, or v normally)
  
  if(!keep.all){
    u1 <- u1[y==1]
    u2 <- u2[y==1]
    d <- d[y==1]
  }
  return(list(N=N, sigma=sigma, B=B, u1=u1, u2=u2, d=d, y=y, N.real=N.real))
}

# obtain a data set by distance sampling a population of N=1000
set.seed(1234)
tmp <-sim.pdata(N=1000, sigma=1, keep.all=FALSE, B=3)
attach(tmp)


# Bin the data and tabulate the bin frequencies. Be sure to pad the 0s!
delta <- 0.5                   # width of distance bins
dist.breaks <-seq(0, B, delta) # make the interval cut points
dclass <- tmp$d %/% delta +1   # Convert distances to categorical distances
nD<-length(dist.breaks)-1     # How many intervals do we have ?
y.obs <- table(dclass)         # next pad the frequency vector
y.padded <- rep(0, nD)
names(y.padded) <- 1:nD
y.padded[names(y.obs)] <- y.obs
y.obs <- y.padded


cp <- c(cp, 1-sum(cp))  # Compute the last cell and add it to the vector

as.vector(rmultinom(n=1, size=1000, prob=cp))


