# TODO: Add comment
# 
# Author: Bruno
###############################################################################



generate_column_names = function(target, r) {
	values = NULL;
	
	for (i in r) {
		values = append( values ,  sprintf ( "%s_%d", target , i) );
	}
	return (values);
}

extract_target = function(fit, target, iterations = NULL){
	
	result = NULL
	
	mainpath = strsplit(target, "_")[[1]]
	## cat("path = ", mainpath, "\n");
	
	if (length(mainpath) == 0) {
		cat("ERROR Invalid variable name(too small)\n",sep = "");
		return (NULL);
	}
	
	variable_name = mainpath[1];
	
	## cat(" - variable_name = ", variable_name, "\n");
	
	if (!is.element(variable_name , names(fit))) {
		cat("ERROR Invalid variable name '",variable_name, "' has not been found)\n",sep = "");
		return (NULL);
	}
	
	# level 1 - get variable
	variable = fit[[variable_name]];
	
	## cat(" - variable taken, class = ",class(variable),"length = ",length(variable),"\n");
	
	explored_iterations = iterations
	if (is.null(explored_iterations)) {
		explored_iterations = c(1:length(variable));
	}
	
	for (iter in explored_iterations) {
		
		path = mainpath[-1];
		
		# level 2 - get iteration
		values = variable[[iter]];
		
		# level 3 - get named_index if list
		
		
		if (is.numeric(values) || is.integer(values)) {
			
			# path can only be a numerical index / we skip
			## cat("   - This is already numerical, we don't need to go down\n")
			
		} else  if (is.list(values))  {
			
			# path is a list we can access an element
			## cat("   - This is a list:", names(values), "\n")
			
			if (length(path) == 0) {
				cat("ERROR Invalid variable name(too small)\n",sep = "");
				return (NULL);
			} 
			
			subitem = path[1];
			path = path[-1];
			
			if (!is.element(subitem , names(values))) {
				cat("ERROR: Cannot find subitem '",subitem,"' among [", paste(names(values), collapse=" "),"]\n",sep = "")
				return (NULL);
			} else {
				values = values[[subitem]];
			}
			
		} else {
			cat("ERROR: UNSUPPORTED TYPE " ,class(values) ," \n",sep = "");
			return (NULL);
		}
		
		
		# now we have values a numerical value or array. 
		# we flatten it or we get what the index requires.
		if (length(path) == 0) {
			## cat("   - We flatten\n");
			# we flatten 
		} else {
			# get the one element selected
			## cat("   - We get the index\n");
			value_index = strtoi(path[1]);
			## cat("   - value_index = " , value_index ,  "\n");
			values = values[value_index];
		}
		
		if (is.null(result)) {
			result = data.frame(t(values))
			names(result) <- generate_column_names(target,c(1:length(values)))
		} else {
			result = rbind(result, values);
		}
	}
	return (result)
}



#'  Extract values within a AM_mcmc_output object. 
#' 
#'  Due to the complexity of AntMAN outputs, AM_mcmc_output object can be difficult
#'  to handle. The AM_extract function ease access of particular variables within 
#'  AM_mcmc_output object.
#'  
#'@param object a \code{\link{AM_mcmc_output}} object
#'@param targets List of variables to extract (ie. K, M, Mna, Tau_mu).
#'@param iterations Can specify particular iterations to extracts, NULL for all.
#'  
#'@export
AM_extract = function(object, targets, iterations = NULL){
	
	df = NULL;
	for (target in targets) {
		
		
		if (target == "CI") {
			## CI Extractor
			nrows = length(object$CI)
			ncols = length(object$CI[[1]])
			tmp = data.frame(array(as.numeric(unlist(object$CI)), dim=c(nrows, ncols)))
			names(tmp) <- generate_column_names(target,c(1:length(ncols)))
		} else {
			## Generic extractor (SLOW)
			tmp = extract_target(object,target,iterations);
		}
		
		
		if (is.null(tmp)) {
			cat("ERROR: Invalid extraction target: ",target,", please make sure this was part of the outputs list in AM_mcmc_parameters.\n", sep="");
			return (NULL);
		}
		
		if(is.null(df)) {df = tmp;}
		else {df = data.frame(df,tmp);}
	}
	return (df);
}

