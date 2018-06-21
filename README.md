# AHMnimble

Run WinBUGS/OpenBUGS/JAGS examples from Volume I of "Applied Hierarchical
Modeling in Ecology" (Kéry and Royle) in NIMBLE.

[*Applied Hierarchical Modeling in Ecology: Analysis of distribution, abundance and species richness in R and BUGS.  Volume I: Prelude and Static Models*](https://www.elsevier.com/books/applied-hierarchical-modeling-in-ecology-analysis-of-distribution-abundance-and-species-richness-in-r-and-bugs/kery/978-0-12-801378-6)
by Marc K&#233;ry and J. Andrew Royle (2015, Academic Press) is a
popular introduction to some common hierarchical models in ecology.  The [book's web site](https://www.mbr-pwrc.usgs.gov/pubanalysis/keryroylebook/) provides code for its examples.

We have adapted most of the examples that use the popular WinBUGS,
OpenBUGS or JAGS packages for Markov chain Monte Carlo (MCMC) to
NIMBLE.  This was easy to do because NIMBLE uses (nearly) the same
modeling language.  We
thank the authors along with Mike Meridith for permission to reproduce
their code with modifications for use with NIMBLE. 

NIMBLE also provides tools for comparing the performance of MCMC
algorithms.   We have generally provided such comparisons between JAGS
and NIMBLE's default sampler configuration.  These were run on a
Mac with a 3.4 GHz Intel iCore 5 with OS X 10.13.4.

We have attempted to make each file labeled with
"*\_example\_nimble.R" or "*\_custom\_nimble.R" be self-contained,
sourcing from other files as needed. 

Some examples remain to be converted.  Please do so and contribute to
this repository.

## Potential to customize NIMBLE
In many cases, it is possible to improve NIMBLE's performance by
re-writing the model using NIMBLE's extensibility and/or customizing
the configuration of MCMC samplers.  As of the initial launch of this
repository, we have only included a few customizations in the MCMC
comparisons.  Please try your own and file a pull request to include
them, or contact us at nimble.stats@gmail.com.

Examples of customizing NIMBLE for Hidden Markov Models (HMMs)
embedded in multi-state capture-recapture models are given in
[Turek et al.](https://link.springer.com/article/10.1007/s10651-016-0353-z).

## How do we compare MCMCs?

We compare MCMCs based on *MCMC efficiency*, which we define as

- MCMC efficiency = Effective sample size / computation time

This carries the interpretation of the number of effectively
independent samples generated per second.

MCMC efficiency is different for each parameter because the effective
sample size is different for each parameter.  For a single metric of
MCMC performance for an entire model, we use the **minimum MCMC
efficiency** across all monitored parameters.  The minimum is
important because one needs to be sure that all parameters have mixed
well before trusting MCMC results. 

Note that we generally do not thin samples when comparing MCMC
efficiency.  One could debate this choice, but we make it because
thinning always results in some loss of statistical information.  In
practice, one often thins.  When comparing MCMC methods, some mix more
slowly but do so computationally faster, and vice-versa.  We attempt a
clean comparison of MCMC efficiency by not thinning.  For some models,
we make an exception and use some thinning to avoid excessively large
MCMC sample sizes.

## Please contribute

Please file pull requests or email us at nimble.stats@gmail.com if you
explore some of these models and have code and/or results that would
be of interest to the community.
