#######################################################################################
###############
############### AntMAN Package
###############
###############
#######################################################################################


#'  compute the hyperparameters of an Normal-Inverse-Gamma distribution using an empirical Bayes approach.
#'  
#'@param y The data y
#'@param scEmu a positive value (default=1) such that marginally E(mu)=(sample variance)*scEmu
#'@param scEsig2 a positive value (default=3) such that marginally E(sig2)=(sample variance)*scEsig2
#'@param CVsig2 The coefficient of variation of sig2 (default=3)
#'  
#'@return hyperparameters
#'@export
AM_emp_bayes_uninorm = function(y,scEmu=1,scEsig2=3,CVsig2=3){
	n <- length(y)   ### sample size
	bary <- mean(y)  ### sample mean
	s2y <- var(y)    ### sample variance
	
	Emu <- bary
	Vmu <- s2y*scEmu
	Esig2 <- s2y/scEsig2
	Vsig2 <- CVsig2^2*Esig2^2
	
	m0    = Emu
	nu0   = 2*(Esig2)^2/Vsig2+4
	sig02 = Esig2*(nu0-2)/nu0
	k0    = sig02/Vmu * nu0/(nu0-2)
	
	return(AM_mix_hyperparams_uninorm(m0=m0,nu0=nu0,sig02=sig02,k0=k0))
	
}

