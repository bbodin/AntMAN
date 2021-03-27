#' Compute the Binder loss function
#' 
#' Compute the Binder loss function of a partitioning
#' with respect to an expected adjacency matrix.
#' 
#' @param eam Expected Adjacency Matrix, i.e., 
#' the matrix whose entries \eqn{E_{ij}} is the posterior probability 
#' that items \eqn{i} and \eqn{j} are together. 
#' If the partitioning is already known, this is just the adjacency matrix.
#' @param labels vector of partition labels. Must be integers. 
#' @param Const_Binder Relative penalty in the Binder loss function 
#' for false-positives vis-a-vis false-negatives. 
#' Must be a real number in the interval [0, 1]. 
#' @return The value of the Binder loss function of the given partition labels
#' with respect to the given pairwise allocation matrix.  
computeBinderLoss <- function(eam, labels, Const_Binder = 0.5){
    if (!is.validAdjacencyMatrix(eam)){
        stop(paste("eam must be a symmetric matrix with values between 0 and 1, ",
                   "and 1s on the diagonal."))
    }
    if (!all(is.finiteInteger(labels))) {
        stop("labels must be a vector of integers.")
    }
    if (length(labels) != ncol(eam)) {
        stop("Incompatible size of eam and labels.")
    }
    return (.computeBinderLoss(eam, labels, Const_Binder))
}

#' Convert partition labels to adjacency matrix
#' 
#' @param labels The integer vector of partition labels.
#' @return A binary integer matrix \eqn{E} where 
#' \eqn{E_{ij} = 1} if \eqn{i} and \eqn{j} if
#' item i and item j have the same partition label, and 0 otherwise.
computeAdjacencyMatrix <- function(labels) {
    if (!all(is.finiteInteger(labels))) {
        stop("labels must be a vector of integers.")
    }
    N <- length(labels)
    return (matrix(as.integer(outer(labels, labels, "==")), ncol=N))
}

#' Binder distance of two partitions
#' 
#' Compute the binder loss function of a partitioning
#' with respect to the adjacency matrix of another partitioning.
#' @param testLabels The vector of integer labels of the partitioning
#' @param refLabels The vector of partition labels to use to construct the adjacency matrix
#' @param Const_Binder Relative penalty in the Binder loss function 
#' for false-positives vis-a-vis false-negatives. 
#' Must be a real number in the interval [0, 1]. 
#' @return The value of the Binder loss function of \code{testLabels} 
#' with respect to the adjacency matrix of \code{refLabels}.
computeBinderDistance <- function(testLabels, refLabels, Const_Binder) {
    if (!all(is.finiteInteger(testLabels))) {
        stop("labels must be a vector of integers.")
    }
    if (!all(is.finiteInteger(refLabels))) {
        stop("labels must be a vector of integers.")
    }
    if (!is.nonNegNumberLessThan1(Const_Binder)) {
        stop("Const_Binder must be a number between 0 and 1.")
    }
    return(computeBinderLoss(computeAdjacencyMatrix(refLabels), testLabels, Const_Binder))
}

is.validAdjacencyMatrix <- function(p) {
    return(all(is.nonNegNumberLessThan1(p)) & 
             nrow(p) == ncol(p) & 
             all(p == t(p)) & 
             sum(diag(p)) == nrow(p))
}

is.finiteInteger <- function(x) {
    return (is.integer(x) & is.finite(x))
}

is.nonNegNumberLessThan1 <- function(x) {
    return (is.numeric(x) & is.finite(x) & x <= 1 & x >= 0)
}