# TODO: Add comment
# 
# Author: Bruno
###############################################################################

library("AntMAN")

weights_prior = AM_mix_weights_prior_gamma(init=2, a=1, b=1)
summary(weights_prior)

mixture_uvn_params = AM_mix_hyperparams_uninorm  (m0=0,k0=0.1,nu0=1,sig02=1.5)
summary(mixture_uvn_params)

mcmc_params        = AM_mcmc_parameters(niter=20000, burnin=1000, thin=10, verbose=1, output = c("ALL"))
summary(mcmc_params)

components_prior   = AM_mix_components_prior_pois (init=3,  a=1, b=1) 
summary(components_prior)
