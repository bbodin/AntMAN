#######################################################################################
###############
############### AntMAN Package
###############
###############
#######################################################################################


#'  Return co-clustering
#'  
#'  Given a MCMC output, this function return co-clustering matrix
#'  
#'@param fit a \code{\link{AM_mcmc_output}} object
#'@return co-clustering matrix
#'
#'@export
AM_coclustering = function (fit) {
	
	result = AM_binder(fit , with_coclustering_probability=TRUE)
	
	return (result[["coclustering_probability"]])
}

