# AntMAN: Anthology of Mixture ANalysis tools

[![Travis build status](https://travis-ci.org/bbodin/AntMAN.svg?branch=master)](https://travis-ci.org/bbodin/AntMAN)

 AntMan is a R package to fit Finite Bayesian Mixture model with random number of component. The MCMC algorithm beyond batman is based on point processes and offer a more computational efficeint alternative to Reversible Jump. Different mixture kernels can be specified: Univariate Gaussian, Univariate Poisson, Univariate Binomial, Multivariate Gaussian, Multivariate Bernoulli (Latent Class Analysis). For the parameters characterising the mixture kernel, we specify conjugate priors, with possibly user specified hyper-parameters. We allow for different choices for the prior on the number of components: Shifted Poisson, Negative Binomial, and Point Masses (i.e. mixtures with fixed number of components).

## How to Install 

### With CRAN 

```
install.packages("AntMAN")
```

### With Github 

First make sure the package is not already installed.

```
remove.packages("AntMAN")
```

Then you can use the following:

```
install.packages("devtools")
devtools::install_github("bbodin/AntMAN", subdir="AntMAN")
```




## How to test

```
library("AntMAN")
data(galaxy, package = "AntMAN")
y_uvn = galaxy
mixture_uvn_params = AntMAN::AM_mix_hyperparams_uninorm  (m0=20.83146, k0=0.3333333, nu0=4.222222, sig02=3.661027)

mcmc_params        = AntMAN::AM_mcmc_parameters(niter=20000, burnin=5000, thin=10, verbose=1)
components_prior   = AntMAN::AM_mix_components_prior_pois (init=3,  a=1, b=1) 
weights_prior      = AntMAN::AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit <- AntMAN::AM_mcmc_fit(
  y = y_uvn,
  mix_kernel_hyperparams = mixture_uvn_params,
  mix_components_prior =components_prior,
  mix_weight_prior = weights_prior,
  mcmc_parameters = mcmc_params)

summary (fit)
plot (fit)
```

## Package Philosophy



The main function of the AntMAN package is ```AM_mcmc_fit``` that performs a Gibbs sampling in order to estimate a mixture of a predifined type ```mix_kernel_hyperparams``` and that models a particular population given a sample ```y```.
Additionnaly a particular component prior ```mix_components_prior``` and a weight prior ```mix_weight_prior``` can be specified and ```mcmc_parameters``` will define the MCMC parameters for the Gibbs sampler (number of interation, burn-in, ...).
A prior on the number of cluster (```init_K```) or a specific allocation (```init_clustering```) can also be specify. Otherwise, the default allocation we assign a different cluster for each element of the sample ```y```. 

For example, in order to identify clusters over a population of patients given a set of medical assumption:

```
mcmc = AM_mcmc_parameters(niter=20000)
mix  = AM_mix_hyperparams_multiber ()
fit = AM_mcmc_fit (mix, mcmc)
summary (fit)
```

In this example ```AM_mix_hyperparams_multiber``` is one of the possible mixture to identify. AntMAN currently support five different mixtures :

```
AM_mix_hyperparams_unipois(alpha0, beta0) 
AM_mix_hyperparams_uninorm(m0, k0, nu0, sig02) 
AM_mix_hyperparams_unibin(a0, b0, mb) 
AM_mix_hyperparams_multiber(a0, b0) 
AM_mix_hyperparams_multinorm(mu0, ka0, nu0, Lam0) 
```

Additionnaly, there is three prior_component available :

```
AM_mix_components_prior_pois
AM_mix_components_prior_negbin
AM_mix_components_prior_dirac
```

For example, in the context of image segmentation, where a maximal number of colour is require, a prior dirac can be used :

```
mcmc = AM_mcmc_parameters(niter=20000)
mix  = AM_mix_hyperparams_multi_norm ()
prior_component = AM_mix_components_prior_dirac(10) # nothing more than 10 colours
fit = AM_mcmc_fit (mix, prior_component, mcmc)
summary (fit)
```

## Acknowledgement 

Thanks for your contributions:
 - David B. Dahl (@dbdahl)
 - Andrea Cremaschi (@AndCre87)

## TODOs

  - Check all the comments in the code
  - File writing feature

## To Remember before submit

  - Please replace \\\% in your Rd-files by \%.
  - Please add and explain the returned objects.
  - Please ensure that your functions do not modify (save or delete) the
user's home filespace in your examples/vignettes/tests. That is not
allow by CRAN policies. Please only write/save files if the user has
specified a directory.
  - Please replace cat() by message() or warning() in your functions (except
for print() and summary() functions). Messages and warnings can be
suppressed if needed.
  - Please always write TRUE and FALSE instead of T and F.
  - Check every headers with wrong authors
  - Please remove the command 'rm(list=ls())' from your files.

