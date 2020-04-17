# TODO: Add comment
# 
# Author: Bruno
###############################################################################




#'  This is a diagnostic function used to compute the log likelihood of MCMC results
#'  
#'  Given a MCMC output, this function return log likelihood of MCMC results
#'  
#'@param MCMC_output_AntMAN a \code{\link{AM_mcmc_output}} object with a full output set (using the ALL argument)
#'  
#'@return log likelihood of MCMC results
#'@export
AM_DIAG_loglikelihood  = function (MCMC_output_AntMAN) {
	
	niter = length(MCMC_output_AntMAN$CI)
	
	loglike_MCMC_iter = rep(0,niter)
	
	for (it in 1:niter) {
		
		loglike_MCMC = 0
		c = MCMC_output_AntMAN$CI[[it]]
		S_m = MCMC_output_AntMAN$S[[it]]
		mu_star = MCMC_output_AntMAN$TAU[[it]]$mu
		Sigma_star = MCMC_output_AntMAN$TAU[[it]]$Sig
		M = MCMC_output_AntMAN$M[[it]]
		U = MCMC_output_AntMAN$U[[it]]
		
		for (i in 1:N) {
			loglike_MCMC = loglike_MCMC + log(S_m[c[i]+1]) + dmvnorm(y[i,], mu_star[c[i]+1,], Sigma_star[,,c[i]+1], TRUE);
		}
		for (j in 1:M) {
			loglike_MCMC = loglike_MCMC + dmvnorm(mu_star[j,], mu0, Sigma_star[,,j]/k0, TRUE);
			loglike_MCMC = loglike_MCMC + log(diwish( Sigma_star[,,j], nu0, Lam0));
			loglike_MCMC = loglike_MCMC + dgamma(S_m[j], gamma_S, scale=1,log= TRUE);
		}
		
		loglike_MCMC = loglike_MCMC + dgamma(Lambda, a2, scale=1/b2, log=TRUE);
		loglike_MCMC = loglike_MCMC + dgamma(gamma_S, a1, scale=1/b1,log= TRUE);
		loglike_MCMC = loglike_MCMC + dgamma(U, M, rate=sum(S_m), log = TRUE);
		
		loglike_MCMC_iter[it] = loglike_MCMC
		
	}
	
	return (loglike_MCMC_iter)
	
}
