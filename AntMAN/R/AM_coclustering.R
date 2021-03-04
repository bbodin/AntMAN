#######################################################################################
###############
############### AntMAN Package
###############
###############
#######################################################################################


#'  Return the co-clustering matrix
#'  
#'
#'  Given an MCMC output produced by \code{AM_MCMC_fit}, this function returns the co-clustering matrix. 
#'  
#'@param fit a \code{\link{AM_mcmc_output}} object
#'@return co-clustering matrix
#'
#'@export
#'
#' @examples
#' ccm <- AM_coclustering(fit)
AM_coclustering = function (fit) {
	
	result = AM_binder(fit , with_coclustering_probability=TRUE)
	
	return (result[["coclustering_probability"]])
}

