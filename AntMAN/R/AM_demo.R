#######################################################################################
###############
############### AntMAN Package
###############
###############
#######################################################################################

#' AM_demo_mvb_prior_poi
#'  
#' Returns AM_MCMC_OUTPUT Object from a running demo of MVB.
#'  
#' @param  n           Number of sample
#' @param  mcmc_params AM_mcmc_parameters object to fit with
#' @keywords demo
#' @export
AM_demo_mvb_poi = function (n = 1000 , mcmc_params = AM_mcmc_parameters(niter=2000, burnin=1000, thin=10, verbose=0, output=c("ALL"))) {
	
	set.seed(123) # TODO : discuss if this make sense
	
	d <- 4
	k <- 3
	TH <- matrix(nrow=k,ncol=d)
	TH[1,] <- c(0.9,0.0,0.2,0.1)
	TH[2,] <- c(0.0,0.9,0.1,0.2)
	TH[3,] <- c(0.0,0.0,0.9,0.9)
	demo_multivariate_binomial <- AM_sample_multibin(n,d,c(0.3,0.3,0.4),TH)
	
	y_mvb  <- demo_multivariate_binomial$y
	ci_mvb <- demo_multivariate_binomial$ci
	
	mixture_mvb_params = AM_mix_hyperparams_multiber  (a0= c(1,1,1,1),b0= c(1,1,1,1))
	
	
	components_prior   = AM_mix_components_prior_pois (init=3,  a=1, b=1) 
	weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)
	
	fit <- AM_mcmc_fit(
			y = y_mvb, 
			mix_kernel_hyperparams = mixture_mvb_params,
			mix_components_prior =components_prior,
			mix_weight_prior = weights_prior,
			mcmc_parameters = mcmc_params)
	
	return (list(input = y_mvb, clusters = ci_mvb, fit = fit))
}


#' AM_demo_mvn_poi
#'  
#' Returns AM_MCMC_OUTPUT Object from a running demo of MVN.
#'  
#' @param  n           Number of sample
#' @param  mcmc_params AM_mcmc_parameters object to fit with
#' @keywords demo
#' @export
AM_demo_mvn_poi = function (n = 1000 , mcmc_params        = AM_mcmc_parameters(niter=4000, burnin=2000, thin=10, verbose=0, output=c("ALL"))) {
	
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
	
	
	
	demo_multivariate_normal <-AM_sample_multinorm(n = n ,d = 2,c(0.3,0.3,0.4),MU,SIG)
	y_mvn  <- demo_multivariate_normal$y
	ci_mvn <- demo_multivariate_normal$ci

	mixture_mvn_params = AM_mix_hyperparams_multinorm   (mu0=c(0,0),ka0=1,nu0=4,Lam0=diag(2))
	
	
	components_prior   = AM_mix_components_prior_pois (init=3,  a=1, b=1) 
	weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)
	
	fit <- AM_mcmc_fit(
			y = y_mvn, 
			mix_kernel_hyperparams = mixture_mvn_params,
			mix_components_prior =components_prior,
			mix_weight_prior = weights_prior,
			mcmc_parameters = mcmc_params)
	
	
	
	return (list(input = y_mvn, clusters = ci_mvn, fit = fit))
}
