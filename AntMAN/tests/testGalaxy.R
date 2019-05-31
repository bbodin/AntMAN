
library("AntMan", lib.loc = "./AntMan.install")

data(galaxy)
y_uvn = galaxy$speed
mixture_uvn_params = AM_uninorm_mix_hyperparams  (m0=0, k0=0.1, nu0=1, sig02=1.5)

mcmc_params        = AM_mcmc_parameters(niter=20000, burnin=5000, thin=10, verbose=1)
components_prior   = AM_mix_components_prior_pois (init=3,  a=1, b=1) 
weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit <- AM_mcmc_fit(
  y = y_uvn,
  mix_kernel_hyperparams = mixture_uvn_params,
  mix_components_prior =components_prior,
  mix_weight_prior = weights_prior,
  mcmc_parameters = mcmc_params)

summary (fit)
plot (fit)
