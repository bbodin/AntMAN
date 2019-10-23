
##############################################
### Load the AntMan package
##############################################

library("AntMAN")

##############################################
### BUILD THE UNIVARIATE BINOMIAL DATA
##############################################

# I fix the random seed
set.seed(123)              
demo_univariate_binomial <-AM_sample_unibin(n=1000,mb=100, pro=c(0.2,0.5,0.3),mth=c(0.1,0.5,0.9)) 
y_uvb  <- demo_univariate_binomial$y
ci_uvb <- demo_univariate_binomial$ci

hist(y_uvb,freq=FALSE,nclass=15,col=colors()[4])
plot(1:length(y_uvb),y_uvb,col=ci_uvb+1)


##############################################################################
### PREPARE THE GIBBS for Binomial mixture with poisson gamma priors
##############################################################################


mixture_uvb_params = AM_unibin_mix_hyperparams  (a0=1,b0=1,N=100)

mcmc_params        = AM_mcmc_parameters(niter=2000, burnin=1000, thin=10, verbose=3)
components_prior   = AM_mix_components_prior_pois (init=3,  a=1, b=1) 
weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit <- AM_mcmc_fit(
       			y = y_uvb, 
                        mix_kernel_hyperparams = mixture_uvb_params,
                        mix_components_prior =components_prior,
                        mix_weight_prior = weights_prior,
                        mcmc_parameters = mcmc_params)

summary (fit)
plot (fit)

