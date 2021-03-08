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
#'  The co-clustering matrix is produced by the simultaneous clustering of the rows and columns. Each entry denotes the (posterior) probability 
#' that items \eqn{i} and \eqn{j} are together.  This technique is also known as
#' bi-clustering and block clustering \insertCite{govaert2013co}{AntMAN}, and is useful for understanding the number of clusters in the dataset. 
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

