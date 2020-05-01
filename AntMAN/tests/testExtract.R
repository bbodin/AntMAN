#######################################################################################
###############
############### AntMAN Package : Tests and Examples
###############
###############
#######################################################################################



if (length(ls())) {
	rm(list = ls());
	.rs.restartR() ;
}


.libPaths( c( "/home/toky/yalenus/research/mixture/AntMan/AntMAN.Rinstall" , .libPaths() ) )
#remove.packages("AntMAN")
library("AntMAN")

n = 100
set.seed(123)


##############################################
### BUILD THE UNIVARIATE NORMAL DATA
##############################################

demo_univariate_normal <-AM_sample_uninorm(n=1000,pro=c(0.2,0.5,0.3),mmu=c(-2.1,0,2.3),ssd=c(0.5,0.5,0.5))
y_uvn  <- demo_univariate_normal$y
ci_uvn <- demo_univariate_normal$ci

hist(y_uvn,freq=FALSE,nclass=15,col=colors()[4])
plot(1:length(y_uvn),y_uvn,col=ci_uvn+1)


##############################################################################
### PREPARE THE GIBBS for Normal mixture with poisson gamma priors
##############################################################################

mixture_uvn_params = AM_mix_hyperparams_uninorm  (m0=0,k0=0.1,nu0=1,sig02=1.5)
mcmc_params        = AM_mcmc_parameters(niter=20000, burnin=1000, thin=10, verbose=0, output = c("ALL"))
components_prior   = AM_mix_components_prior_pois (init=3,  a=1, b=1) 
weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)

uvn_fit <- AM_mcmc_fit(
		y = y_uvn, 
		mix_kernel_hyperparams = mixture_uvn_params,
		mix_components_prior =components_prior,
		mix_weight_prior = weights_prior,
		mcmc_parameters = mcmc_params)

summary (uvn_fit)



##############################################################################
### TEST EXTRACT
##############################################################################


df_uvn = AM_extract(uvn_fit,c("CI","mu","sig2","W","PREDICTIVE","U","M","K","gamma","lambda"))
CI_uvn = AM_extract(uvn_fit,c("CI"))





##############################################
### BUILD THE MULTIVARIATE NORMAL DATA
##############################################

set.seed(123)


MU <- matrix(nrow=3,ncol=2)

MU[1,] <- c(0,0)
MU[2,] <- c(-3,-3)
MU[3,] <- c(4,4)


sig1 <- c(1,1)
rho1 <- 0
Sig1 <- matrix(c(sig1[1]^2,rho1*sig1[1]*sig1[2], rho1*sig1[1]*sig1[2],sig1[2]^2),byrow=TRUE,nrow=2) 

sig2 <- c(1,1)
rho2 <- -0.7
Sig2 <- matrix(c(sig2[1]^2,rho2*sig2[1]*sig2[2], rho2*sig2[1]*sig2[2],sig2[2]^2),byrow=TRUE,nrow=2) 

sig3 <- c(1,1)
rho3 <- -0.3
Sig3 <- matrix(c(sig3[1]^2,rho3*sig3[1]*sig3[2], rho3*sig3[1]*sig3[2],sig3[2]^2),byrow=TRUE,nrow=2) 


SIG <- array(0,dim=c(3,2,2))
SIG[1,,] <- Sig1
SIG[2,,] <- Sig2
SIG[3,,] <- Sig3



demo_multivariate_normal <-AM_sample_multinorm(n = 1000 ,d = 2,c(0.3,0.3,0.4),MU,SIG)
y_mvn  <- demo_multivariate_normal$y
ci_mvn <- demo_multivariate_normal$ci

hist(y_mvn,freq=FALSE,nclass=15,col=colors()[4])
plot(y_mvn,col=ci_mvn+1)


##############################################################################
### PREPARE THE GIBBS for multivariate Normal mixture with poisson gamma priors
##############################################################################


mixture_mvn_params = AM_mix_hyperparams_multinorm   (mu0=c(0,0),ka0=1,nu0=4,Lam0=diag(2))

mcmc_params        = AM_mcmc_parameters(niter=4000, burnin=2000, thin=10, verbose=0, output=c("ALL"))
components_prior   = AM_mix_components_prior_pois (init=3,  a=1, b=1) 
weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit_mvn <- AM_mcmc_fit(
		y = y_mvn, 
		mix_kernel_hyperparams = mixture_mvn_params,
		mix_components_prior =components_prior,
		mix_weight_prior = weights_prior,
		mcmc_parameters = mcmc_params)


summary (fit_mvn)


##############################################################################
### TEST EXTRACT
##############################################################################



CI_mvn = AM_extract(fit_mvn,c("CI"))
df_mvn = AM_extract(fit_mvn,c("mu","Sig","W","U","M","K","gamma","lambda"))





##############################################
### BUILD THE UNIVARIATE POISSON DATA
##############################################

set.seed(123)
demo_univariate_poisson <-AM_sample_unipois(n=1000,pro=c(0.2,0.5,0.3),mth=c(5,25,50)) 
y_uvp  <- demo_univariate_poisson$y
ci_uvp <- demo_univariate_poisson$ci

hist(y_uvp,freq=FALSE,nclass=15,col=colors()[4])
plot(1:length(y_uvp),y_uvp,col=ci_uvp+1)



##############################################################################
### PREPARE THE GIBBS for Poisson mixture with poisson gamma priors
##############################################################################

fit_poisson <- AM_mcmc_fit(
		y = y_uvp, initial_clustering = 0:(length(y_uvp)-1),
		mix_kernel_hyperparams = AM_mix_hyperparams_unipois (alpha0=2, beta0=0.2),
		mix_components_prior =AM_mix_components_prior_pois (init=3,  a=1, b=1) ,
		mix_weight_prior =AM_mix_weights_prior_gamma(init=2, a=1, b=1),
		mcmc_parameters = AM_mcmc_parameters(niter=2000, burnin=1000, thin=10, verbose=0, output=c("ALL"))
)



##############################################################################
### TEST EXTRACT
##############################################################################

CI_poisson = AM_extract(fit_poisson,c("CI"))
df_poisson = AM_extract(fit_poisson,c("theta","W","U","M","K","gamma","lambda"))


