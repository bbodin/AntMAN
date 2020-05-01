#######################################################################################
###############
############### AntMAN Package
###############
###############
#######################################################################################


#' plot AM_mcmc_output  
#' 
#' plot some useful informations about the mcmc results
#'  
#'@param x a AM_mcmc_output object
#'@param tags A list of variables to consider
#'@param title Title for the plot
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
		if (!is.null(x$H) && (length(x$H) > 0) && (length(x$H[[1]]) > 0)) {targets = c(targets,paste0("H_",names(x$H[[1]])))}
		if (!is.null(x$Q) && (length(x$Q) > 0) && (length(x$Q[[1]]) > 0)) {targets = c(targets,paste0("Q_",names(x$Q[[1]])))}
	}
	message("Plotting pairs from ",paste(targets,collapse=","));
	
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
#'@param title Title for the plot
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
		if (!is.null(x$H) && (length(x$H) > 0) && (length(x$H[[1]]) > 0)) {targets = c(targets,paste0("H_",names(x$H[[1]])))}
		if (!is.null(x$Q) && (length(x$Q) > 0) && (length(x$Q[[1]]) > 0)) {targets = c(targets,paste0("Q_",names(x$Q[[1]])))}
	}
	message("Plotting density from ",paste(targets,collapse=","));
	
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
#'@param title Title for the plot
#'@return TBD
#'@importFrom bayesplot mcmc_trace color_scheme_set
#'@export
AM_plot_traces=function(x,tags = NULL,title = "MCMC Results"){
	
	targets = tags
	if (is.null(targets)) {
		if (!is.null(x$M)) {targets = c(targets,"M")}
		if (!is.null(x$M)) {targets = c(targets,"K")}
		if (!is.null(x$MNA)) {targets = c(targets,"MNA")}
		if (!is.null(x$H) && (length(x$H) > 0) && (length(x$H[[1]]) > 0)) {targets = c(targets,paste0("H_",names(x$H[[1]])))}
		if (!is.null(x$Q) && (length(x$Q) > 0) && (length(x$Q[[1]]) > 0)) {targets = c(targets,paste0("Q_",names(x$Q[[1]])))}
	}
	message("Plotting traces from ",paste(targets,collapse=","));
	
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
#'@param title Title for the plot
#'@return TBD
#'@importFrom bayesplot mcmc_intervals color_scheme_set
#'@export
AM_plot_values=function(x,tags = NULL,title = "MCMC Results"){
	
	targets = tags
	if (is.null(targets)) {
		if (!is.null(x$M)) {targets = c(targets,"M")}
		if (!is.null(x$M)) {targets = c(targets,"K")}
		if (!is.null(x$MNA)) {targets = c(targets,"MNA")}
		if (!is.null(x$H) && (length(x$H) > 0) && (length(x$H[[1]]) > 0)) {targets = c(targets,paste0("H_",names(x$H[[1]])))}
		if (!is.null(x$Q) && (length(x$Q) > 0) && (length(x$Q[[1]]) > 0)) {targets = c(targets,paste0("Q_",names(x$Q[[1]])))}
	}
	message("Plotting values from ",paste(targets,collapse=","));
	
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
	
	sorted = TRUE
	arguments <- list(...)
	if ("sorted" %in% names(arguments)) {
		sorted = arguments[[sorted]]
	}
	
	
	if (!is.null(x$CI)) {
		
		message("Plotting Similarity Matrix");
		binder_result = AM_binder(x , with_coclustering_probability=TRUE)
		
		clustering = binder_result[["clustering"]]
		pij = binder_result[["coclustering_probability"]]
		
		if (sorted) {
			new_indexes = order(clustering, decreasing = T)
			pij = pij[new_indexes,new_indexes];
		}
		
		image(pij,main="Similarity matrix",col = gray.colors(30))
	} else {
		warning("CI has not been generated. Cannot plot the similarity matrix.")
	}
	
}
#'  Plot the Autocorrelation
#'  
#'  Given a MCMC output, this function will produce acf bar
#'  
#'@param x a AM_mcmc_output object
#'@param tags A list of variables to consider
#'@param title Title for the plot
#'  
#'@importFrom bayesplot mcmc_acf_bar
#'@export
AM_plot_chaincor=function(x, tags = NULL, title = "MCMC Results"){
	
	targets = tags
	if (is.null(targets)) {
		if (!is.null(x$M)) {targets = c(targets,"M")}
		if (!is.null(x$M)) {targets = c(targets,"K")}
		if (!is.null(x$MNA)) {targets = c(targets,"MNA")}
		if (!is.null(x$H) && (length(x$H) > 0) && (length(x$H[[1]]) > 0)) {targets = c(targets,paste0("H_",names(x$H[[1]])))}
		if (!is.null(x$Q) && (length(x$Q) > 0) && (length(x$Q[[1]]) > 0)) {targets = c(targets,paste0("Q_",names(x$Q[[1]])))}
	}
	message("Plotting Autocorrelation from ",paste(targets,collapse=","))
	
	if (length(targets) > 0) {
		df = AM_extract(x,targets)
		color_scheme_set("brightblue")
		mcmc_acf_bar(df)
	}
	

}
	
