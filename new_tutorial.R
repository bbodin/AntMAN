
##############################################
### Load the AntMan package
##############################################
library("AntMan", lib.loc = "./AntMan.install")


##############################################
### BUILD THE UNIVARIATE POISSON DATA
##############################################
rm(list=ls())
set.seed(123)
demo_univariate_poisson <-AM_sample_unipois(n=1000,pro=c(0.2,0.5,0.3),mth=c(5,25,50)) 
y_uvp  <- demo_univariate_poisson$y
ci_uvp <- demo_univariate_poisson$ci

hist(y_uvp,freq=F,nclass=15,col=colors()[4])
plot(1:length(y_uvp),y_uvp,col=ci_uvp+1)



##############################################################################
### PREPARE THE GIBBS for Poisson mixture with poisson gamma priors
##############################################################################
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




##############################################
### BUILD THE UNIVARIATE NORMAL DATA
##############################################

rm(list=ls())
set.seed(123)
demo_univariate_normal <-AM_sample_uninorm(n=1000,pro=c(0.2,0.5,0.3),mmu=c(-2.1,0,2.3),ssd=c(0.5,0.5,0.5))
y_uvn  <- demo_univariate_normal$y
ci_uvn <- demo_univariate_normal$ci

hist(y_uvn,freq=F,nclass=15,col=colors()[4])
plot(1:length(y_uvn),y_uvn,col=ci_uvn+1)



##############################################################################
### PREPARE THE GIBBS for Normal mixture with dirac gamma priors
##############################################################################

mixture_uvn_params = AM_uninorm_mix_hyperparams  (m0=0,k0=0.1,nu0=1,sig02=1.5)

mcmc_params        = AM_mcmc_parameters(niter=2000, burnin=1000, thin=10, verbose=1)
components_prior   = AM_mix_components_prior_dirac (Mstar=3) 
weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit <- AM_mcmc_fit(
       			y = y_uvn, init_K=1,
                        mix_kernel_hyperparams = mixture_uvn_params,
                        mix_components_prior =components_prior,
                        mix_weight_prior = weights_prior,
                        mcmc_parameters = mcmc_params)

summary (fit)
plot (fit)


##############################################
### BUILD THE UNIVARIATE NORMAL DATA
##############################################

rm(list=ls())
set.seed(123)
demo_univariate_normal <-AM_sample_uninorm(n=1000,pro=c(0.2,0.5,0.3),mmu=c(-2.1,0,2.3),ssd=c(0.5,0.5,0.5))
y_uvn  <- demo_univariate_normal$y
ci_uvn <- demo_univariate_normal$ci

hist(y_uvn,freq=F,nclass=15,col=colors()[4])
plot(1:length(y_uvn),y_uvn,col=ci_uvn+1)



mixture_uvn_params = AM_uninorm_mix_hyperparams  (m0=0,k0=0.1,nu0=1,sig02=1.5)

mcmc_params        = AM_mcmc_parameters(niter=2000, burnin=1000, thin=10, verbose=1)
components_prior   = AM_mix_components_prior_negbin (M_P=0.1, M_R=1) 
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





##############################################
### BUILD THE UNIVARIATE POISSON DATA
##############################################
rm(list=ls())
set.seed(123)
demo_univariate_poisson <-AM_sample_unipois(n=1000,pro=c(0.2,0.5,0.3),mth=c(5,25,50)) 
y_uvp  <- demo_univariate_poisson$y
ci_uvp <- demo_univariate_poisson$ci

hist(y_uvp,freq=F,nclass=15,col=colors()[4])
plot(1:length(y_uvp),y_uvp,col=ci_uvp+1)



##############################################################################
### [NEW INTERFACE] PREPARE THE GIBBS for Poisson mixture with poisson dirac priors
##############################################################################
mcmc_params        = AM_mcmc_parameters(niter=2000, burnin=1000, thin=10, verbose=3)
mixture_uvp_params = AM_unipois_mix_hyperparams (alpha0=2, beta0=0.2)
components_prior   = AM_mix_components_prior_dirac (Mstar=5) 
weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)
init_ci_uvp <- 0:(length(y_uvp)-1);

fit_poisson_dirac <- AM_mcmc_fit(
       			y = y_uvp, init_K=1,
                        mix_kernel_hyperparams = mixture_uvp_params,
                        mix_components_prior =components_prior,
                        mix_weight_prior = weights_prior,
                        mcmc_parameters = mcmc_params)


summary (fit_poisson_dirac)
plot (fit_poisson_dirac)

##############################################
### BUILD THE UNIVARIATE NORMAL DATA
##############################################

rm(list=ls())
set.seed(123)
demo_univariate_normal <-AM_sample_uninorm(n=1000,pro=c(0.2,0.5,0.3),mmu=c(-2.1,0,2.3),ssd=c(0.5,0.5,0.5))
y_uvn  <- demo_univariate_normal$y
ci_uvn <- demo_univariate_normal$ci

hist(y_uvn,freq=F,nclass=15,col=colors()[4])
plot(1:length(y_uvn),y_uvn,col=ci_uvn+1)



##############################################################################
### PREPARE THE GIBBS for Normal mixture with poisson gamma priors
##############################################################################

mixture_uvn_params = AM_uninorm_mix_hyperparams  (m0=0,k0=0.1,nu0=1,sig02=1.5)

mcmc_params        = AM_mcmc_parameters(niter=2000, burnin=1000, thin=10, verbose=3)
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



##############################################
### BUILD THE UNIVARIATE BINOMIAL DATA
##############################################

# I fix the random seed
rm(list=ls())
set.seed(123)              
demo_univariate_binomial <-AM_sample_unibin(n=1000,mb=100, pro=c(0.2,0.5,0.3),mth=c(0.1,0.5,0.9)) 
y_uvb  <- demo_univariate_binomial$y
ci_uvb <- demo_univariate_binomial$ci

hist(y_uvb,freq=F,nclass=15,col=colors()[4])
plot(1:length(y_uvb),y_uvb,col=ci_uvb+1)


##############################################################################
### PREPARE THE GIBBS for Binomial mixture with poisson gamma priors
##############################################################################


mixture_uvb_params = AM_unibin_mix_hyperparams  (a0=1,b0=1,mb=100)

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


##############################################
### BUILD THE MULTIVARIATE NORMAL DATA
##############################################

rm(list=ls())
set.seed(123)


MU <- matrix(nrow=3,ncol=2)

MU[1,] <- c(0,0)
MU[2,] <- c(-3,-3)
MU[3,] <- c(4,4)


sig1 <- c(1,1)
rho1 <- 0
Sig1 <- matrix(c(sig1[1]^2,rho1*sig1[1]*sig1[2], rho1*sig1[1]*sig1[2],sig1[2]^2),byrow=T,nrow=2) 

sig2 <- c(1,1)
rho2 <- -0.7
Sig2 <- matrix(c(sig2[1]^2,rho2*sig2[1]*sig2[2], rho2*sig2[1]*sig2[2],sig2[2]^2),byrow=T,nrow=2) 

sig3 <- c(1,1)
rho3 <- -0.3
Sig3 <- matrix(c(sig3[1]^2,rho3*sig3[1]*sig3[2], rho3*sig3[1]*sig3[2],sig3[2]^2),byrow=T,nrow=2) 


SIG <- array(0,dim=c(3,2,2))
SIG[1,,] <- Sig1
SIG[2,,] <- Sig2
SIG[3,,] <- Sig3


library(mvtnorm) ##  TODO Why do I need this !!!!
demo_multivariate_normal <-AM_sample_multinorm(n = 1000 ,d = 2,c(0.3,0.3,0.4),MU,SIG)
y_mvn  <- demo_multivariate_normal$y
ci_mvn <- demo_multivariate_normal$ci

hist(y_mvn,freq=F,nclass=15,col=colors()[4])
plot(y_mvn,col=ci_mvn+1)


##############################################################################
### PREPARE THE GIBBS for multivariate Normal mixture with poisson gamma priors
##############################################################################


mixture_mvn_params = AM_multinorm_mix_hyperparams   (mu0=c(0,0),ka0=1,nu0=4,Lam0=diag(2))

mcmc_params        = AM_mcmc_parameters(niter=2000, burnin=1000, thin=10, verbose=3)
components_prior   = AM_mix_components_prior_pois (init=3,  a=1, b=1) 
weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit <- AM_mcmc_fit(
       			y = y_mvn, 
                        mix_kernel_hyperparams = mixture_mvn_params,
                        mix_components_prior =components_prior,
                        mix_weight_prior = weights_prior,
                        mcmc_parameters = mcmc_params)

summary (fit)
plot (fit)

##############################################
### BUILD THE MULTIVARIATE BINOMIAL DATA
##############################################

rm(list=ls())
set.seed(123)


d <- 4
k <- 3
TH <- matrix(nrow=k,ncol=d)

TH[1,] <- c(0.9,0.0,0.2,0.1)
TH[2,] <- c(0.0,0.9,0.1,0.2)
TH[3,] <- c(0.0,0.0,0.9,0.9)
n <- 1000
demo_multivariate_binomial <- AM_sample_multibin(n,d,c(0.3,0.3,0.4),TH)

y_mvb  <- demo_multivariate_binomial$y
ci_mvb <- demo_multivariate_binomial$ci

hist(y_mvb,freq=F,nclass=15,col=colors()[4])
plot(y_mvb,col=ci_mvb+1)


##############################################################################
### PREPARE THE GIBBS for multivariate BINOMIAL mixture with poisson gamma priors
##############################################################################


mixture_mvb_params = AM_multiber_mix_hyperparams  (a0= c(1,1,1,1),b0= c(1,1,1,1))

mcmc_params        = AM_mcmc_parameters(niter=2000, burnin=1000, thin=10, verbose=3)
components_prior   = AM_mix_components_prior_pois (init=3,  a=1, b=1) 
weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit <- AM_mcmc_fit(
       			y = y_mvb, 
                        mix_kernel_hyperparams = mixture_mvb_params,
                        mix_components_prior =components_prior,
                        mix_weight_prior = weights_prior,
                        mcmc_parameters = mcmc_params)

summary (fit)
plot (fit)

1