# TODO: Add comment
# 
# Author: Bruno
###############################################################################



#'  Run the binder algorithm using R (TBD)
#'  
#'  TBD TODO : need the C version of that.
#'  
#'@param fit                            Output from MCMC_fit
#'@param weight                         Weight between bad and good pairs, default is 0.5.
#'@param with_coclustering_probability  By default this function only return the index of the closest guess. 
#'                                      When with_coclustering_probability, the function also return the coclustering probability matrix.
#'  
#'@export
AM_binder = function (fit,  weight = 0.5, with_coclustering_probability = FALSE) {
	CI = data.matrix(AM_extract(fit, c("CI")));
	#Equal costs
	Const_Binder <- weight
	
	#c_out contains matrix of allocation labels
	c_out <- CI
	
	
	#c_contains matrix of allocation labels
	
	
	n_save <- dim(c_out)[1]
	N      <- dim(c_out)[2]
	
	
	
	#Compute similarity matrix pij
	
	pij <- matrix(0,N,N)
	
	for(g in 1:n_save){
		
		pij <- pij + outer(c_out[g,], c_out[g,], "==")
		
	}
	
	pij <- pij/n_save
	
	
	
	Binder_f <- rep(0,n_save)
	
	for(g in 1:n_save){
		
		cij <- outer(c_out[g,], c_out[g,], "==")
		
		aux <- (pij - Const_Binder) * as.matrix(cij)
		
		aux <-  aux[upper.tri(aux)]
		
		Binder_f[g] <- sum(aux)
		
	}
	
	Binder_ind <- which.max(Binder_f)
	
	if (with_coclustering_probability) {
		return (list(coclustering_probability = pij, clustering = CI[Binder_ind,], index = Binder_ind));
	} else {
		return (list(clustering = CI[Binder_ind,], index = Binder_ind));
	}
}