# Notes on examples from AHM chapter 5.

This chapter focuses on simple linear, generalized linear and generalized linear mixed models (LM, GLM, and GLMM).

For the simplest models, JAGS is more efficient than NIMBLE.  As the chapter progresses, NIMBLE is more efficient than JAGS on more complicated models.  All of the models mix easily with both packages until Section 5.14, where we see a slow-mixing GLMM.  NIMBLE outperforms JAGS on that example.
