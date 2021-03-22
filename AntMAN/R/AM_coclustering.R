#######################################################################################
###############
############### AntMAN Package
###############
###############
#######################################################################################


#'  Return the co-clustering matrix
#'  
#'
#'  Given an MCMC output produced by \link{\code{AM_MCMC_fit}}, this function returns the co-clustering matrix. 
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
	

	CI =  as.matrix((AM_extract(fit, c("CI")))[["CI"]]);

	n_save <- dim(CI)[1]
	N  <- dim(CI)[2]
	
	c_out <- lapply(1:n_save, function(x){outer(CI[x,], CI[x,], "==")})
	pij = Reduce("+", c_out)/n_save

	#cast pij into matrix
	pij = matrix(unlist(pij), ncol=N, nrow=N, byrow=F)


	# result = AM_binder(fit , with_coclustering_probability=TRUE)
	
	# return (result[["coclustering_probability"]])
	
	# #c_out contains matrix of allocation labels
	# c_out <- CI
	
	
	# #c_contains matrix of allocation labels
	
	
	# n_save <- dim(c_out)[1]
	# N      <- dim(c_out)[2]
	
	
	# #Compute similarity matrix pij
	
	# pij <- matrix(0,N,N)
	
	# for(g in 1:n_save){
		
	# 	pij <- pij + outer(c_out[g,], c_out[g,], "==")
		
	# }
	
	# pij <- pij/n_save
	
	return(coclustering_probability = pij);
	#return (list(coclustering_probability = pij));
}

