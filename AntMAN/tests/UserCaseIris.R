#######################################################################################
###############
############### AntMAN Package : Tests and Examples
###############
###############
#######################################################################################

##############################################
### Load the AntMan package
##############################################

library("AntMAN")
set.seed(123)

##############################################
### Load the IRIS package
##############################################


data("iris")

mat <- matrix(0, 150, 4)
mat[,1] <- iris$Sepal.Length
mat[,2] <- iris$Sepal.Width
mat[,3] <- iris$Petal.Length
mat[,4] <- iris$Petal.Width


mixture_mvn_params <- AM_mix_hyperparams_multinorm(mu0 = rep(0, 4),
                                                   ka0 = 1,
                                                   nu0 = 6,
                                                   Lam0= 1*diag(4))
mcmc_params<-  AM_mcmc_parameters(niter=20000,
                                  burnin=5000,
                                  thin=10,
                                  verbose=1)


components_prior <- AM_mix_components_prior_pois (a=10,b=2)
weights_prior <- AM_mix_weights_prior_gamma(a=1, b=1)


fit <- AM_mcmc_fit(
  y = mat,
  mix_kernel_hyperparams = mixture_mvn_params,
  mix_components_prior =components_prior,
  mix_weight_prior = weights_prior,
  mcmc_parameters = mcmc_params)
