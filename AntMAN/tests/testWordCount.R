quit() ## Skip this test - too long
library("AntMAN")

data(said)

y_uvp = said$Frequency
  
mcmc_params        = AM_mcmc_parameters(niter=2000, burnin=1000, thin=10, verbose=3, output=c("CI","K"))
mixture_uvp_params = AM_unipois_mix_hyperparams (alpha0=2, beta0=0.2)
components_prior   = AM_mix_components_prior_pois (init=3,  a=1, b=1) 
weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)
init_ci_uvp <- 0:(length(y_uvp)-1);

fit_poisson <- AM_mcmc_fit(
			y = y_uvp, initial_clustering = init_ci_uvp,
                        mix_kernel_hyperparams = mixture_uvp_params,
                        mix_components_prior =components_prior,
                        mix_weight_prior = weights_prior,
                        mcmc_parameters = mcmc_params)


summary (fit_poisson)
plot (fit_poisson)
