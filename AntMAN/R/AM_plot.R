# TODO: Add comment
# 
# Author: Bruno
###############################################################################


#' plot AM_mcmc_output  
#' 
#' plot some useful informations about the mcmc results
#'  
#'@param x a AM_mcmc_output object
#'@param tags A list of variables to pairplot
#'@return Same as ggpairs function, a ggmatrix object that if called, will print.
#'@importFrom graphics image
#'@importFrom GGally ggpairs
#'@importFrom grDevices gray.colors
#'@export
AM_plot_pairs=function(x,tags = NULL,title = "MCMC Results"){
	
	targets = tags
	if (is.null(targets)) {
		if (!is.null(x$M)) {targets = c(targets,"M")}
		if (!is.null(x$M)) {targets = c(targets,"K")}
		if (!is.null(x$MNA)) {targets = c(targets,"MNA")}
		if (!is.null(x$H)) {targets = c(targets,paste0("H_",names(x$H[[1]])))}
		if (!is.null(x$Q)) {targets = c(targets,paste0("Q_",names(x$Q[[1]])))}
		#if (!is.null(x$S)) {targets = c(targets,"S")}
		#if (!is.null(x$U)) {targets = c(targets,"U")}
		#if (!is.null(x$TAU)) {targets = c(targets,paste0("TAU_",names(x$TAU[[1]])))}
	}
	message("Will plot informations from ",paste(targets,collapse=","))
	
	if (length(targets) > 0) {
		df = AM_extract(x,targets)
		ggpairs(df, title = title, upper = list(continuous = "points"),)
	}

}


#' plot density of variables from AM_mcmc_output object  
#' 
#' TBD
#'  
#'@param x a AM_mcmc_output object
#'@param tags A list of variables to consider
#'@return TBD
#'@importFrom bayesplot mcmc_areas color_scheme_set
#'@export
AM_plot_density=function(x,tags = NULL,title = "MCMC Results"){
	## TODO / Need raffa script for integer values
	## skip integer
	targets = tags
	if (is.null(targets)) {
		if (!is.null(x$M)) {targets = c(targets,"M")}
		if (!is.null(x$M)) {targets = c(targets,"K")}
		if (!is.null(x$MNA)) {targets = c(targets,"MNA")}
		if (!is.null(x$H)) {targets = c(targets,paste0("H_",names(x$H[[1]])))}
		if (!is.null(x$Q)) {targets = c(targets,paste0("Q_",names(x$Q[[1]])))}
		#if (!is.null(x$S)) {targets = c(targets,"S")}
		#if (!is.null(x$U)) {targets = c(targets,"U")}
		#if (!is.null(x$TAU)) {targets = c(targets,paste0("TAU_",names(x$TAU[[1]])))}
	}
	message("Will plot informations from ",paste(targets,collapse=","))
	
	if (length(targets) > 0) {
		df = AM_extract(x,targets)
		color_scheme_set("brightblue")
		mcmc_areas(df)
	}
	
}

#' plot traces of variables from AM_mcmc_output object  
#' 
#' TBD
#'  
#'@param x a AM_mcmc_output object
#'@param tags A list of variables to consider
#'@return TBD
#'@importFrom bayesplot mcmc_trace color_scheme_set
#'@export
AM_plot_traces=function(x,tags = NULL,title = "MCMC Results"){
	
	targets = tags
	if (is.null(targets)) {
		if (!is.null(x$M)) {targets = c(targets,"M")}
		if (!is.null(x$M)) {targets = c(targets,"K")}
		if (!is.null(x$MNA)) {targets = c(targets,"MNA")}
		if (!is.null(x$H)) {targets = c(targets,paste0("H_",names(x$H[[1]])))}
		if (!is.null(x$Q)) {targets = c(targets,paste0("Q_",names(x$Q[[1]])))}
		#if (!is.null(x$S)) {targets = c(targets,"S")}
		#if (!is.null(x$U)) {targets = c(targets,"U")}
		#if (!is.null(x$TAU)) {targets = c(targets,paste0("TAU_",names(x$TAU[[1]])))}
	}
	message("Will plot informations from ",paste(targets,collapse=","))
	
	if (length(targets) > 0) {
		df = AM_extract(x,targets)
		color_scheme_set("brightblue")
		mcmc_trace(df)
	}
	
}

#' plot histogram of variables from AM_mcmc_output object  
#' 
#' TBD
#'  
#'@param x a AM_mcmc_output object
#'@param tags A list of variables to consider
#'@return TBD
#'@importFrom bayesplot mcmc_intervals color_scheme_set
#'@export
AM_plot_values=function(x,tags = NULL,title = "MCMC Results"){
	
	targets = tags
	if (is.null(targets)) {
		if (!is.null(x$M)) {targets = c(targets,"M")}
		if (!is.null(x$M)) {targets = c(targets,"K")}
		if (!is.null(x$MNA)) {targets = c(targets,"MNA")}
		if (!is.null(x$H)) {targets = c(targets,paste0("H_",names(x$H[[1]])))}
		if (!is.null(x$Q)) {targets = c(targets,paste0("Q_",names(x$Q[[1]])))}
		#if (!is.null(x$S)) {targets = c(targets,"S")}
		#if (!is.null(x$U)) {targets = c(targets,"U")}
		#if (!is.null(x$TAU)) {targets = c(targets,paste0("TAU_",names(x$TAU[[1]])))}
	}
	message("Will plot informations from ",paste(targets,collapse=","))
	
	if (length(targets) > 0) {
		df = AM_extract(x,targets)
		color_scheme_set("brightblue")
		mcmc_intervals(df)
	}
	
}


#'  Plot the Similarity Matrix
#'  
#'  Given a MCMC output, this function will produce an image of the Similarity Matrix
#'  
#'@param x a AM_mcmc_output object
#'@param ... All additional parameters wil lbe pass to the image command.
#'  
#'@export
AM_plot_similarity_matrix=function(x, ...){
	if (!is.null(x$CI)) {
		res = AM_coclustering(x)
		image(res,main="Similarity matrix",col = gray.colors(30))
	} else {
		warning("CI has not been generated. Cannot plot the similarity matrix.")
	}
	
}
