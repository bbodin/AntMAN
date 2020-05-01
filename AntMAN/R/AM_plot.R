#######################################################################################
###############
############### AntMAN Package
###############
###############
#######################################################################################

density_discrete_variables <- function(Par, color=rgb(0.4, 0.8, 1, alpha=0.7), single_maxy=TRUE, ...){
	rows <- dim(Par)[2]
	fun <- function(xx){
		return(table(xx)/length(xx))
	}
	
	tables <- lapply(Par,fun)
	#cat("maxy is",maxy,"\n")
	if(single_maxy){
		maxy=zeroes(rows)
		for(r in 1:rows){
			maxy[r]=as.numeric(max(tables[[r]]))
		}
		#	print("ciao")
	}else{
		maxy <- rep(max(unlist(tables)),rows)
		#	print("CIAO")
	}
	#print(maxy)
	### Allow the user to set a 
	
	par(mfrow=c(rows,1))
	for(r in 1:rows){
		plot(tables[[r]],lwd=8,col=color,ylim=c(0,maxy[r]),xlab=names(tables)[r],ylab="p.m.f.",...)
	}
	
	
	
}


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
		if (!is.null(x$K)) {targets = c(targets,"K")}
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
		if (!is.null(x$K)) {targets = c(targets,"K")}
	}
	message("Plotting density from ",paste(targets,collapse=","));
	
	if (length(targets) > 0) {
		df = AM_extract(x,targets)
		color_scheme_set("brightblue")
		mcmc_areas(df)
	}
	
}


#' TBD
#' 
#' TBD
#'  
#'@param x a AM_mcmc_output object
#'@param tags A list of variables to consider
#'@param title Title for the plot
#'@return TBD
#'@importFrom grDevices rgb
#'@export
AM_plot_bars=function(x,tags = NULL,title = "MCMC Results"){
	## TODO / Need raffa script for integer values
	## skip integer
	targets = tags
	if (is.null(targets)) {
		if (!is.null(x$M)) {targets = c(targets,"M")}
		if (!is.null(x$K)) {targets = c(targets,"K")}
	}
	message("Plotting density from ",paste(targets,collapse=","));
	
	if (length(targets) > 0) {
		df = AM_extract(x,targets)
		density_discrete_variables(df, single_maxy=TRUE)
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
		if (!is.null(x$K)) {targets = c(targets,"K")}
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
		if (!is.null(x$K)) {targets = c(targets,"K")}
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
	
