#######################################################################################
###############
############### AntMAN Package
###############
###############
#######################################################################################

#################################################################################
##### Package Definition
#################################################################################

#' AntMAN: A package for fitting Finite Bayesian Mixture model with random number of component.
#'
#'@description AntMAN: Anthology of Mixture ANalysis tools 
#' AntMan is a R package to fit Finite Bayesian Mixture  model with random number of component. 
#' The MCMC algorithm beyond AntMan is based on point processes and offer a more computational  
#' efficient  alternative to Reversible Jump. 
#' Different mixture kernels can be specified: Univariate Gaussian, Univariate Poisson, Univariate Binomial, Multivariate Gaussian, 
#' Multivariate Bernoulli (Latent Class Analysis). For the parameters characterising the mixture kernel, we specify 
#' conjugate priors, with possibly user specified hyper-parameters.   
#' We allow for different choices for the prior on the number of components: 
#' Shifted Poisson, Negative Binomial, and Point Masses (i.e. mixtures with fixed number of components).
#' 
#'@section Prior functions:
#' The Prior functions ...
#' 
#'@section Package Philosophy:
#' 
#' The main function of the AntMAN package is \code{\link{AM_mcmc_fit}}. AntMAN  performs a Gibbs sampling in order to fit, 
#' in a Bayesian framework, a mixture model of a predifined type \code{mix_kernel_hyperparams}  given a sample \code{y}. 
#' Additionally AntMAN allows the user to specify a prior on the number of components \code{mix_components_prior} and on the weights  \code{mix_weight_prior} of the mixture.
#' MCMC parameters \code{mcmc_parameters} need to be given as argument for the Gibbs sampler (number of interation, burn-in, ...). 
#' Initial values for the number of cluster (\code{init_K}) or a specific clustering allocation (\code{init_clustering}) can also be user-specify. 
#' Otherwise, by the default allocation we assign a different cluster for each element of the sample \code{y} as initial allocation. This choice can be computetionally inefficient. 
#' 
#' 
#' For example, in order to identify clusters over a population of patients given a set of medical assumptions:
#' 
#'```
#' mcmc = AM_mcmc_parameters(niter=20000) 
#' mix = AM_mix_hyperparams_multiber () 
#' fit = AM_mcmc_fit (mix, mcmc) 
#' summary (fit)
#'```
#' 
#' In this example \code{AM_mix_hyperparams_multiber} is one of the possible mixture to identify. 
#' 
#' AntMAN currently support five different mixtures :
#' 
#' ```
#' AM_mix_hyperparams_unipois(alpha0, beta0) 
#' AM_mix_hyperparams_uninorm(m0, k0, nu0, sig02) 
#' AM_mix_hyperparams_unibin(a0, b0, mb) 
#' AM_mix_hyperparams_multiber(a0, b0) 
#' AM_mix_hyperparams_multinorm(mu0, ka0, nu0, Lam0)
#' ```
#' 
#' Additionnaly, there is three prior_component available :
#' 
#' ```
#' AM_mix_components_prior_pois()
#' AM_mix_components_prior_negbin() 
#' AM_mix_components_prior_dirac()
#' ```
#' 
#' For example, in the context of image segmentation, where a maximal number of colour is require, a prior dirac can be used :
#' 
#' ```
#' mcmc = AM_mcmc_parameters(niter=20000) 
#' mix = AM_mix_hyperparams_multinorn () 
#' prior_component = AM_mix_components_prior_dirac(10) # nothing more than 10 colours 
#' fit = AM_mcmc_fit (mix, prior_component, mcmc) summary (fit)
#' ```
#' 
#'@importFrom Rcpp evalCpp
#'@importFrom stats kmeans rbinom rnorm rpois runif sd acf density quantile var
#'@importFrom graphics plot hist rasterImage abline layout legend lines
#'@importFrom salso dlso
#'@docType package
#'@name AntMAN
NULL




