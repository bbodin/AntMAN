# TODO: Add comment
# 
# Author: Bruno
###############################################################################



#'  Run the binder algorithm using R (TBD)
#'  
#'  TBD
#'  
#'@param C fit CI
#'@param Const_Binder Const_Binder
#'  
#'@export
AM_binder=function (CI,  weight = 0.5) {
	
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
	
	
	
	return (Binder_ind)
}