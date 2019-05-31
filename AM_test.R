##############################################
### Load the AntMan package
##############################################
library("AntMAN", lib.loc = "./AntMAN.Rinstall")


set.seed(123)
test_priors();

rm(list=ls())
set.seed(123)
demo_univariate_normal <-AM_sample_uninorm(n=1000,pro=c(0.2,0.5,0.3),mmu=c(-2.1,0,2.3),ssd=c(0.5,0.5,0.5))
y_uvn  <- demo_univariate_normal$y
ci_uvn <- demo_univariate_normal$ci

hist(y_uvn,freq=F,nclass=15,col=colors()[4])
plot(1:length(y_uvn),y_uvn,col=ci_uvn+1)



mixture_uvn_params = AM_uninorm_mix_hyperparams  (m0=0,k0=0.1,nu0=1,sig02=1.5)

mcmc_params        = AM_mcmc_parameters(niter=2000, burnin=1000, thin=10, verbose=1)
components_prior   = AM_mix_components_prior_negbin (M_P=0.9, M_R=1) 
weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit <- AM_mcmc_fit(
       			y = y_uvn, init_K=1,
                        mix_kernel_hyperparams = mixture_uvn_params,
                        mix_components_prior =components_prior,
                        mix_weight_prior = weights_prior,
                        mcmc_parameters = mcmc_params)

summary (fit)
plot(fit$K_post)
plot(y_uvn,col=fit$ci_post[[length(fit$ci_post)]]+1)
