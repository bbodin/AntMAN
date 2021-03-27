#######################################################################################
###############
############### AntMAN Package
###############
###############
#######################################################################################


#' Sequentially Allocated Latent Structure Optimisation
#'
#' Heuristic partitioning to minimise the expected Binder loss function 
#' with respect to a given expected adjacency matrix. 
#'
#' @param fit an \code{\link{AM_mcmc_fit}} object. 
#' @param maxClusts Maximum number of clusters. 
#' The actual number of clusters searched may be lower. 
#' If set to 0L, the maximum is automatically limited by the number of items. 
#' @param Const_Binder Relative penalty in the Binder loss function 
#' for false-positives vis-a-vis false-negatives. 
#' Must be a real number in the interval [0, 1]. 
#' @param batchSize Number of permutations scanned per thread. 
#' If set to 0L, the thread will continue to scan permutations until it times out 
#' (in which case \code{timeLimit} cannot be 0L).
#' @param nScans Number of scans for each permutation. 
#' @param maxThreads Maximum number of threads to use. 
#' If set to 0L (default), the maximum number of threads 
#' will be determined by the runtime. 
#' Set to 1L for no parallelisation. 
#' The actual number of threads used may be lower than \code{maxThreads}.
#' @param timeLimit Maximum computation time for each thread in milliseconds. 
#' The actual computational time may be higher, 
#' since the time limit is only checked at the end of each iteration. 
#' If set to 0L, the thread will never time out 
#' (in which case \code{batchSize} cannot be 0L).
#' @return A list containing the following items:
#' \itemize{
#' \item \code{Labels} - the vector of partition labels
#' \item \code{BinderLoss} - the associated binder loss function
#' \item \code{NumClusts} - the number of clusters found
#' \item \code{NumPermutations} - the number of permutations actually scanned
#' \item \code{WallClockTime} - cumulative wall-clock time used by all threads in milliseconds
#' \item \code{TimeLimitReached} - whether the computation time limit was reached in any of the threads
#' \item \code{NumThreads} - actual number of threads used.
#' }
#'@export
AM_binder <- function(fit, maxClusts=0L, Const_Binder = 0.5, batchSize = 1000L, nScans = 10L, maxThreads = 0L, timeLimit = 600000L) {

	eam  = AM_coclustering(fit)
		
	if (!is.validAdjacencyMatrix(eam)){
		stop(paste("eam must be a symmetric matrix with values between 0 and 1, ",
						"and 1s on the diagonal."))
	}
	if (!is.finiteInteger(maxClusts) | 
			maxClusts < 0) {
		stop("maxClusts must be a nonnegative integer. ")
	}
	if (!is.nonNegNumberLessThan1(Const_Binder)) {
		stop("Const_Binder must be a number between 0 and 1.")
	}
	if (!is.finiteInteger(batchSize) | 
			batchSize < 0) {
		stop("batchSize must be a nonnegative integer.")
	}
	if (!is.finiteInteger(nScans) | 
			nScans <= 0) {
		stop("nScans must be a positive integer.")
	}
	if (!is.finiteInteger(maxThreads) | 
			maxThreads < 0) {
		stop("maxThreads must be a nonnegative integer.")
	}
	if (!is.finiteInteger(timeLimit) | 
			timeLimit < 0) {
		stop("timeLimit must be a nonnegative integer.")
	}
	if (timeLimit == 0L & batchSize == 0L) {
		stop("batchSize and timeLimit cannot both be 0.")
	}
	
	temp = .salso(eam, maxClusts, Const_Binder, batchSize, nScans, maxThreads, timeLimit)
	temp[["Labels"]] <- as.integer(as.vector(temp[["Labels"]]))
	return (temp)
}


