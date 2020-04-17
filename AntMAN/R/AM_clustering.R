#######################################################################################
###############
############### AntMAN Package
###############
###############
#######################################################################################

#'  Return maximum likelihood estimation (laugreen)
#'  
#'  Given a MCMC output, this function return maximum likelihood estimation.
#'  
#'@param fit a \code{\link{AM_mcmc_output}} object
#'@param C used to speed up the function, can be the coclustering as returned from the command \code{\link{AM_coclustering}}.
#'@return maximum likelihood estimation (laugreen)
#'@export
AM_clustering_estimation_laugreen = function (fit, C = NULL) {
	FF <- vector("numeric")
	K <- 0.5
	G <- length(fit$K)
	n = length(fit$CI[[1]])
	if (C == NULL) C = AM_coclustering(fit) 
	ci <- t(do.call(cbind,fit$CI))+1
	for(g in 1:(G)){
		ss <- ci[g,]
		cij <- outer(ss,ss,'==')
		pluto <- (C-K)*as.matrix(cij)
		pluto <-  pluto[upper.tri(pluto)]
		FF[g] <- sum(pluto)
	}
	ind.bind <- which.max(FF)[1]
	
	clust_bind <- ci[ind.bind,]
	return (clust_bind)
}


#'  Return maximum likelihood estimation (squared_loss)
#'  
#'  Given a MCMC output, this function return maximum likelihood estimation.
#'  
#'@param fit a \code{\link{AM_mcmc_output}} object
#'  
#'@importFrom salso dlso psm
#'@return maximum likelihood estimation (squared_loss)
#'@export
AM_clustering_estimation_squared_loss = function (fit) {
	if(requireNamespace("salso")) {
		draws = t(do.call(cbind,fit$CI))
		fres = dlso(psm(draws),"binder",draws)
		return (fres);
	} else {
		stop ( "The package salso is required when using this command." );
	}
}

#'  Return maximum likelihood estimation (average)
#'  
#'  Given a MCMC output, this function return maximum likelihood estimation.
#'  
#'@param fit a \code{\link{AM_mcmc_output}} object
#'@return maximum likelihood estimation (average)
#'  
#'@importFrom mcclust minbinder comp.psm
#'@export
AM_clustering_estimation_average = function (fit) {
	if(requireNamespace("mcclust")) {
		mcinput = t(do.call(cbind,fit$CI))
		psm2 <- comp.psm(mcinput+1)
		mbind2 <- minbinder(psm2)
		names(mbind2)
		return (mbind2$cl)
	} else {
		stop ( "The package mcclust is required when using this command." );
	}
}




#'  Return co-clustering
#'  
#'  Given a MCMC output, this function return co-clustering matrix
#'  
#'@param fit a \code{\link{AM_mcmc_output}} object
#'@return co-clustering matrix
#'
#'@export
AM_coclustering = function (fit) {
	
	G <- length(fit$K)
	n = length(fit$CI[[1]])
	C <- matrix(0,ncol=n,nrow=n)
	ci <- t(do.call(cbind,fit$CI))+1
	for(g in 1:G){
		ss <- ci[g,]
		cij <- outer(ss,ss,'==')
		C <- C + cij
	}
	
	return ( C / G )
}

#'  Return co-clustering slowly
#'  
#'  Given a MCMC output, this function return co-clustering matrix
#'  
#'@param fit a \code{\link{AM_mcmc_output}} object
#'@return co-clustering matrix
#'  
#'
#'@export
AM_coclustering_slow = function (fit) {
	
	n = length(fit$CI[[1]])
	res = matrix(0,n,n)
	for (i in (fit$CI)) {
		drow = matrix(rep(i,each=n),nrow=n)
		dcol = matrix(rep(i,each=n), ncol=n, byrow=TRUE)	
		res = res + !(drow-dcol)
	}
	return (res / max(res))
}

