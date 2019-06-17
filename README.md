# antman

What has to be done : 
  - DESCRIPTION file 
  - Check all the comments in the code
  - File writing feature
  - Correct Parallel arguments
  


The main function of the AntMAN package is \code{AM_mcmc_fit} that performs a Gibbs sampling in order to estimate a mixture of a predifined type \code{mix_kernel_hyperparams} and that models a particular population given a sample \code{y}.
Additionnaly a particular component prior \code{mix_components_prior} and a weight prior \code{mix_weight_prior} can be specified and \code{mcmc_parameters} will define the MCMC parameters for the Gibbs sampler (number of interation, burn-in, ...).
A prior on the number of cluster (\code{init_K}) or a specific allocation (\code{init_clustering}) can also be specify. Otherwise, the default allocation we assign a different cluster for each element of the sample \code{y}. 

For example, in order to identify clusters over a population of patients given a set of medical assumption:

mcmc = AM_mcmc_parameters(niter=20000)
mix  = AM_multiber_mix_hyperparams ()
fit = AM_mcmc_fit (mix, mcmc)
summary (fit)

In this example \code{AM_multiber_mix_hyperparams} is one of the possible mixture to identify. AntMAN currently support five different mixtures :

AM_unipois_mix_hyperparams(alpha0, beta0) 
AM_uninorm_mix_hyperparams(m0, k0, nu0, sig02) 
AM_unibin_mix_hyperparams(a0, b0, mb) 
AM_multiber_mix_hyperparams(a0, b0) 
AM_multinorm_mix_hyperparams(mu0, ka0, nu0, Lam0) 

Additionnaly, there is three prior_component available :

AM_mix_components_prior_pois
AM_mix_components_prior_negbin
AM_mix_components_prior_dirac

For example, in the context of image segmentation, where a maximal number of colour is require, a prior dirac can be used :

mcmc = AM_mcmc_parameters(niter=20000)
mix  = AM_multinorn_mix_hyperparams ()
prior_component = AM_mix_components_prior_dirac(10) # nothing more than 10 colours
fit = AM_mcmc_fit (mix, prior_component, mcmc)
summary (fit)