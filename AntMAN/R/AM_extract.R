#######################################################################################
###############
############### AntMAN Package
###############
###############
#######################################################################################


generate_column_names = function(target, r) {
	values = NULL;
	
	for (i in r) {
		values = append( values ,  sprintf ( "%s_%d", target , i) );
	}
	return (values);
}



extract_target = function(fit, target, iterations = NULL, debug = FALSE){
	
	result = NULL
	
	mainpath = strsplit(target, "_")[[1]]
	if (debug) message("path = ", paste(mainpath,collapse = " "), "\n");
	
	if (length(mainpath) == 0) {
		warning("ERROR Invalid variable name(too small)\n");
		return (NULL);
	}
	
	variable_name = mainpath[1];
	
	if (debug) message(" - variable_name = ", variable_name, "\n");
	
	if (!is.element(variable_name , names(fit))) {
		warning("ERROR Invalid variable name '",variable_name, "' has not been found)\n");
		return (NULL);
	}
	
	# level 1 - get variable
	variable = fit[[variable_name]];
	
	if (debug) message(" - variable taken, class = ",class(variable),", length = ",length(variable),"\n");
	
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
			if (debug) message("   - This is already numerical, we don't need to go down\n")
			
		} else  if (is.list(values))  {
			
			# path is a list we can access an element
			if (debug) message("   - This is a list:", paste(names(values), collapse = " "), "\n")
			
			if (length(path) == 0) {
				warning("ERROR Invalid variable name(too small)\n");
				return (NULL);
			} 
			
			subitem = path[1];
			path = path[-1];
			
			if (!is.element(subitem , names(values))) {
				warning("ERROR: Cannot find subitem '",subitem,"' among [", paste(names(values), collapse=" "),"]\n")
				return (NULL);
			} else {
				if (debug) message("   - We access subitem: ", subitem, "\n")
				values = values[[subitem]];
			}
			
		} else {
			warning("ERROR: UNSUPPORTED TYPE " ,class(values) ," \n");
			return (NULL);
		}
		
		
		if (is.numeric(values) || is.integer(values)) {
			if (debug) message("   - lower level can be handle\n");
		} else {
			warning("ERROR: UNSUPPORTED TYPE " ,class(values) ," \n");
			return (NULL);
		}
		
		# now we have values a numerical value or array. 
		# we flatten it or we get what the index requires.
		if (length(path) == 0) {
			if (debug) message("   - We need to flatten\n");
			# we flatten 
		} else {
			# get the one element selected
			if (debug) message("   - We get the index\n");
			value_index = strtoi(path[1]);
			if (debug) message("   - value_index = " , value_index ,  "\n");
			values = values[value_index];
		}
		
		values = as.vector(unlist(values))
		
		if (is.null(result)) {
			result = data.frame(t(values))
			names(result) <- generate_column_names(target,c(1:length(values)))
		} else {
			if (ncol(result) > length(values)) {
				values = c( values , rep(NA, ncol(result) - length(values))); 
			}
			
			while (ncol(result) < length(values)) { ## TODO : Please find more efficient!
				result = cbind(result,c(NA))
				names(result) <- generate_column_names(target,c(1:length(values)))
			}
			if (ncol(result) != length(values)) {
				warning("ERROR: NUMBER of COLUMN CHANGED FROM " ,ncol(result) , " to ", length(values)," \n");
				return (NULL);
			} else {
				result = rbind(result, values);
			}
		}
	}
	return (result)
}