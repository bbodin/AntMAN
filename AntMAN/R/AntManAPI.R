
#################################################################################
##### Package Definition
#################################################################################

#' AntMAN: A package for fitting Finite Bayesian Mixture model with random number of component.
#'
#'@description AntMAN: Anthology of Mixture ANalysis tools 
#' AntMan is a R package to fit Finite Bayesian Mixture  model with random number of component. 
#' The MCMC algorithm beyond AntMan is based on point processes and offer a more computational  
#' efficeint  alternative to Reversible Jump. 
#' Different mixture kernels can be specified: Univariate Gaussian, Univariate Poisson, Univariate Binomial, Multivariate Gaussian, 
#' Multivariate Bernoulli (Latent Class Analysis). For the parameters characterising the mixture kernel, we specify 
#' conjugate priors, with possibly user specified hyper-parameters.   
#' We allow for different choices for the prior on the number of components: 
#' Shifted Poisson, Negative Binomial, and Point Masses (i.e. mixtures with fixed number of components).
#' 
#'@section Prior functions:
#' The Prior functions ...
#' 
#'@section Package Philosophy:
#' 
#' The main function of the AntMAN package is \code{\link{AM_mcmc_fit}}. AntMAN  performs a Gibbs sampling in order to fit, 
#' in a Bayesian framework, a mixture model of a predifined type \code{mix_kernel_hyperparams}  given a sample \code{y}. 
#' Additionally AntMAN allows the user to specify a prior on the number of components \code{mix_components_prior} and on the weights  \code{mix_weight_prior} of the mixture.
#' MCMC parameters \code{mcmc_parameters} need to be given as argument for the Gibbs sampler (number of interation, burn-in, ...). 
#' Initial values for the number of cluster (\code{init_K}) or a specific clustering allocation (\code{init_clustering}) can also be user-specify. 
#' Otherwise, by the default allocation we assign a different cluster for each element of the sample \code{y} as initial allocation. This choice can be computetionally inefficient. 
#' 
#' 
#' For example, in order to identify clusters over a population of patients given a set of medical assumptions:
#' 
#'```
#' mcmc = AM_mcmc_parameters(niter=20000) 
#' mix = AM_multiber_mix_hyperparams () 
#' fit = AM_mcmc_fit (mix, mcmc) 
#' summary (fit)
#'```
#' 
#' In this example \code{AM_multiber_mix_hyperparams} is one of the possible mixture to identify. 
#' 
#' AntMAN currently support five different mixtures :
#' 
#' ```
#' AM_unipois_mix_hyperparams(alpha0, beta0) 
#' AM_uninorm_mix_hyperparams(m0, k0, nu0, sig02) 
#' AM_unibin_mix_hyperparams(a0, b0, mb) 
#' AM_multiber_mix_hyperparams(a0, b0) 
#' AM_multinorm_mix_hyperparams(mu0, ka0, nu0, Lam0)
#' ```
#' 
#' Additionnaly, there is three prior_component available :
#' 
#' ```
#' AM_mix_components_prior_pois()
#' AM_mix_components_prior_negbin() 
#' AM_mix_components_prior_dirac()
#' ```
#' 
#' For example, in the context of image segmentation, where a maximal number of colour is require, a prior dirac can be used :
#' 
#' ```
#' mcmc = AM_mcmc_parameters(niter=20000) 
#' mix = AM_multinorn_mix_hyperparams () 
#' prior_component = AM_mix_components_prior_dirac(10) # nothing more than 10 colours 
#' fit = AM_mcmc_fit (mix, prior_component, mcmc) summary (fit)
#' ```
#' 
#'@importFrom Rcpp evalCpp
#'@importFrom stats kmeans rbinom rnorm rpois runif sd acf density quantile var
#'@importFrom graphics plot hist rasterImage abline layout legend lines
#'@importFrom sdols dlso
#'@importFrom corrplot corrplot
#'@docType package
#'@name AntMAN
NULL


#################################################################################
##### Datasets Definition
#################################################################################


#' Usage frequency of the word said in the Brown corpus.
#'     
#'@format  A list with 500 observations on the frequency of said in different texts.
#'  
#'@source https://www.kaggle.com/nltkdata/brown-corpus
#'@references  Francis, W., and Kucera, H. (1982) Frequency Analysis of English Usage, Houghton Mifflin Company, Boston.
#'@examples
#'  data(said)
#'@keywords datasets
#'@docType data
"said"

#' Galaxy velocities dataset
#' 
#' This data set consider physical information on velocities (km/second) for
#' 82 galaxies reported by Roeder (1990). These are drawn from six well-separated
#' conic sections of the Corona Borealis region.
#'  
#'@format A data frame with X rows and Y variables.
#'  
#'@format A numeric vector giving the speed of galaxies (1000*(km/second)) 
#'@examples
#'  data(galaxy)
#'@source  Roeder, K. (1990) Density estimation with confidence sets exemplified by superclusters and voids in the galaxies, Journal of the American Statistical Association, 85: 617-624.  
#'@keywords datasets
#'@docType data
"galaxy"


#' carcinoma
#'  
#' The carcinoma data from Agresti (2002, 542) consist of seven dichotomous variables that represent 
#' the ratings by seven pathologists of 118 slides on the presence or absence of carcinoma.
#' The purpose of studying these data is to model "interobserver agreement" by examining how
#' subjects might be divided into groups depending upon the consistency of their diagnoses.
#'  
#'@format A data frame with 118 rows and 7 variables (from A to G).
#'  
#'@references Agresti A (2002). Categorical Data Analysis. John Wiley & Sons, Hoboken.
#'  
#'@examples
#'  data(carcinoma)
#'  
#'@keywords datasets
#'@docType data
"carcinoma"

#' Teen Brain Images from the National Institutes of Health, U.S.
#'
#' Picture of brain activities from a teenager consumming drugs. 
#'  
#'@format  A list that contains \code{dim} a (W:width,H:height) pair, and  \code{pic} a data frame (W*H pixels image in RGB format).
#'  
#'@source https://www.flickr.com/photos/nida-nih/29741916012
#'@references Crowley TJ, Dalwani MS, Mikulich-Gilbertson SK, Young SE, Sakai JT, Raymond KM, et al. (2015) Adolescents' Neural Processing of Risky Decisions: Effects of Sex and Behavioral Disinhibition. PLoS ONE 10(7): e0132322. doi:10.1371/journal.pone.0132322
#'  
#'@examples
#'  data(brain)
#'  
#'@keywords datasets
#'@docType data
"brain"


#################################################################################
##### AM_mcmc_output, Summary, Plot.
#################################################################################


#'  Plot the Similarity Matrix
#'  
#'  Given a MCMC output, this function will produce an image of the Similarity Matrix
#'  
#'@param fit a \code{\link{AM_mcmc_output}} object
#'  
#'@export
AM_plot_similarity_matrix=function(C, ...){
	#### Co clustering probability : How many time two are in the same custer. 

	n = dim(C)[1]
	## library(corrplot)
	## col3 <- colorRampPalette(c("red", "white", "blue")) 
	## corrplot(res, diag = FALSE, method = "color", type = "upper",col = col3(100), cl.lim = c(0, 1), tl.pos = "n")
	image(1:n,1:n,C,xaxt='n', yaxt="n",main="Similarity matrix",col = gray.colors(30), ...)
}

#' S3 class AM_mcmc_output.
#' @description Output type of return values from  \code{\link{AM_mcmc_fit}}. See paper for the moment. 
#' @exportClass AM_mcmc_output
#' @seealso \code{\link{AM_mcmc_fit}}
#' @name AM_mcmc_output
NULL


#' plot AM_mcmc_output  
#' 
#' plot some useful informations about the mcmc results
#'  
#'@param x a AM_mcmc_output object
#'  
#'@method plot AM_mcmc_output 
#'@export
plot.AM_mcmc_output=function(x,...){
  if (!is.null(x$M))  {hist(x$M,main="M values") ; readline(prompt="Press [enter] to continue");}
  if (!is.null(x$K))  {hist(x$K,main="K Values") ; readline(prompt="Press [enter] to continue");}
  
  #### Co clustering probability : How many time two are in the same custer. 
  if (!is.null(x$CI)) {
	  G <- length(x$K)
	  n = length(x$CI[[1]])
	  res = AM_coclustering(x)
	  ## library(corrplot)
	  ## col3 <- colorRampPalette(c("red", "white", "blue")) 
	  ## corrplot(res, diag = FALSE, method = "color", type = "upper",col = col3(100), cl.lim = c(0, 1), tl.pos = "n")
	  image(1:n,1:n,res,xaxt='n', yaxt="n",main="Similarity matrix",col = gray.colors(30))
	  readline(prompt="Press [enter] to continue");
  }
}

#'  summary AM_mcmc_output 
#'  
#'  Print some useful informations about the mcmc results
#'  
#'@param object a \code{\link{AM_mcmc_output}} object
#'  
#'  
#'@method summary AM_mcmc_output 
#'@export
summary.AM_mcmc_output=function(object,...){
	cat( "Name", "\t", "Mean","\t", "StdDev","\n");
	if (!is.null(object$K)) cat( "K", "\t", mean(object$K),"\t", sd(object$K),"\n");
	if (!is.null(object$M)) cat( "M", "\t", mean(object$M),"\t", sd(object$M),"\n");
	if (!is.null(object$Mna)) cat( "Mna", "\t", mean(object$Mna),"\t", sd(object$Mna),"\n");
}
#'  Return maximum likelihood estimation (laugreen)
#'  
#'  Given a MCMC output, this function return maximum likelihood estimation.
#'  
#'@param fit a \code{\link{AM_mcmc_output}} object
#'  
#'@export
AM_clustering_estimation_laugreen = function (fit, C = NULL) {
	FF <- vector("numeric")
	K <- 0.5
	G <- length(fit$K)
	n = length(fit$CI[[1]])
	if (C == NULL) C = AM_coclustering(fit) 
	ci <- t(do.call(cbind,fit$CI))+1
	for(g in 1:(G)){
		ss <- ci[g,]
		cij <- outer(ss,ss,'==')
		pluto <- (C-K)*as.matrix(cij)
		pluto <-  pluto[upper.tri(pluto)]
		FF[g] <- sum(pluto)
	}
	ind.bind <- which.max(FF)[1]
	
	clust_bind <- ci[ind.bind,]
	return (clust_bind)
}


#'  Return maximum likelihood estimation (squared_loss)
#'  
#'  Given a MCMC output, this function return maximum likelihood estimation.
#'  
#'@param fit a \code{\link{AM_mcmc_output}} object
#'  
#'@importFrom sdols dlso
#'@export
AM_clustering_estimation_squared_loss = function (fit) {
	library('sdols')
	fres = dlso(t(do.call(cbind,fit$CI)))
	return (fres);
}

#'  Return maximum likelihood estimation (average)
#'  
#'  Given a MCMC output, this function return maximum likelihood estimation.
#'  
#'@param fit a \code{\link{AM_mcmc_output}} object
#'  
#'@importFrom mcclust minbinder comp.psm
#'@export
AM_clustering_estimation_average = function (fit) {
	library("mcclust")
	mcinput = t(do.call(cbind,fit$CI))
	psm2 <- comp.psm(mcinput+1)
	mbind2 <- minbinder(psm2)
	names(mbind2)
	return (mbind2$cl)
}




#'  Return co-clustering
#'  
#'  Given a MCMC output, this function return co-clustering matrix
#'  
#'@param fit a \code{\link{AM_mcmc_output}} object
#'  
#'
#'@export
AM_coclustering = function (fit) {
	
	G <- length(fit$K)
	n = length(fit$CI[[1]])
	C <- matrix(0,ncol=n,nrow=n)
	ci <- t(do.call(cbind,fit$CI))+1
	for(g in 1:G){
		ss <- ci[g,]
		cij <- outer(ss,ss,'==')
		C <- C + cij
	}
	
	return ( C / G )
}

#'  Return co-clustering slowly
#'  
#'  Given a MCMC output, this function return co-clustering matrix
#'  
#'@param fit a \code{\link{AM_mcmc_output}} object
#'  
#'
#'@export
AM_coclustering_slow = function (fit) {

	n = length(fit$CI[[1]])
	res = matrix(0,n,n)
	for (i in (fit$CI)) {
		drow = matrix(rep(i,each=n),nrow=n)
		dcol = matrix(rep(i,each=n), ncol=n, byrow=TRUE)	
		res = res + !(drow-dcol)
	}
	return (res / max(res))
}

#################################################################################
##### AM_mcmc_fit function
#################################################################################


#' Performs a Gibbs sampling
#' 
#' The \code{AM_mcmc_fit} function performs a Gibbs sampling in order to estimate 
#' a mixture of a predifined type \code{mix_kernel_hyperparams} (defined with \code{AM_*_mix_hyperparams} functions, where star 
#' \code{*} denotes the chosen kernel) 
#' sample data \code{y}. 
#' Additionaly a prior distribution on the number of mixture components  
#' must be specified  through  \code{mix_components_prior} 
#' (generated with  \code{AM_mix_components_prior_*} functions, where \code{*}  denotes the chosen prior). Similarly,  
#' a prior on the weights of  the mixture is specified through \code{mix_weight_prior} 
#' (defined with  \code{AM_mix_weights_prior_*} functions). Finally, with \code{mcmc_parameters} the user sets
#' the MCMC parameters for the Gibbs sampler (defined with  \code{AM_mcmc_parameters} functions). 
#' 
#' If no initial clustering is specified (either as \code{init_K} or \code{init_clustering}), 
#' then every observation is allocated to a different cluster. 
#' If \code{init_K} is specified then AntMAN initialises the clustering through K-means. 
#' 
#' **Warning**: if the user does not specify init_K or initial_cluster, the first steps can be be time-consuming because of default value for the initial clustering. 
#' 
#'
#'@param y input data, can be a vector or a matrix.
#'@param mix_kernel_hyperparams is a configuration list, defined by *_mix_hyperparams functions, where * denotes the chosen kernel.
#'@param initial_clustering is a vector CI of initial cluster assignement. If no clustering is specified (either as \code{init_K} or \code{init_clustering}), then every observation is assigned to its own cluster.
#'@param fixed_clustering is a vector CI of cluster assignement that will remained unchanged for every iterations.
#'@param init_K initial value for  the number of cluster. When this is specified, AntMAN intitialises the clustering assign usng K-means.
#'@param mix_components_prior is a configuration list defined by AM_mix_components_prior_* functions, where * denotes the chosen prior.
#'@param mix_weight_prior is a configuration list defined by AM_weight_prior_* functions, where * denotes the chosen prior specification.
#'@param mcmc_parameters is a configuration list defined by AM_mcmc_parameters. 
#'@return The return value is a \code{\link{AM_mcmc_output}} object.
#'@examples
#' AM_mcmc_fit( AM_sample_unipois()$y, 
#'              AM_unipois_mix_hyperparams (alpha0=2, beta0=0.2), 
#'              mcmc_parameters = AM_mcmc_parameters(niter=2000, burnin=1000, thin=10))
#'@useDynLib AntMAN
#'@export

AM_mcmc_fit <- function(
  y, 
  mix_kernel_hyperparams, 
  initial_clustering = NULL, 
  init_K = NULL, 
  fixed_clustering = NULL, 
  mix_components_prior = AM_mix_components_prior_pois() , 
  mix_weight_prior = AM_mix_weights_prior_gamma(), 
  mcmc_parameters = AM_mcmc_parameters() ) {
  fixed_cluster = FALSE
  if (is.null(fixed_clustering) & is.null(init_K) & !is.null(initial_clustering)) {
	  fixed_cluster = FALSE
  } else if (!is.null(init_K) & is.null(initial_clustering)& is.null(fixed_clustering)) {
	  fixed_cluster = FALSE
    initial_clustering <- kmeans(y, init_K)$cluster
  } else if (is.null(init_K) & is.null(initial_clustering)& is.null(fixed_clustering)) {
	  fixed_cluster = FALSE
    initial_clustering <- 0:(length(y)-1)
  } else if (is.null(init_K) & is.null(initial_clustering)& !is.null(fixed_clustering)) { 
	  fixed_cluster = TRUE
	  initial_clustering = fixed_clustering
  } 
  else {
    stop("Please provide only one of K_init or initial_clustering or fixed_clustering.")
  }
  
  structure(IAM_mcmc_fit(y = y, mix_kernel_hyperparams = mix_kernel_hyperparams, initial_clustering = initial_clustering, fixed_clustering=  fixed_cluster , mix_components_prior = mix_components_prior, mix_weight_prior = mix_weight_prior, mcmc_parameters = mcmc_parameters)
            , class = "AM_mcmc_output") 
}

#################################################################################
##### AM_mcmc_parameters function
#################################################################################



#' MCMC Parameters
#' 
#' This function generates an MCMC parameters list to be used as \code{mcmc_parameters} argument within \code{AM_mcmc_fit}. 
#' 
#' The \code{niter} argument specify the total number of iteration. 
#' \code{burnin} is the number of initial iterations to discard.
#' \code{thin} specifies how often a draw from teh posterior distribution is stored after 
#' burnin, i.e. one every -th samples is saved. Therefore, the toral number of MCMC samples saved is 
#' (\code{niter} -\code{burnin})/\code{thin}. If thin =1, then AntMAN stores every iteration.  
#' 
#' 
#'@param niter        Total number of MCMC iterations. 
#'@param burnin       Number of iterations to discard as burn-in.
#'@param thin         Thining rate.
#'@param verbose      A value from 0 to 4, that specifies the desired level of verbosity (0:None, 1:Warnings, 2:Debug, 3:Extras)
#'@param output       A list of parameters output to return
#'@param output_dir   Path to an output dir, when to store all the outputs.
#'@param parallel     Some of the algorithms can be run in parallel using OpenMP. When set to True, this parameter triggers the parallelism.
#'@return list to be used as \code{mcmc_parameters} argument with \code{AM_mcmc_fit}. 
#'@examples 
#' AM_mcmc_parameters (niter=1000, burnin=10000, thin=50)
#' AM_mcmc_parameters (niter=1000, burnin=10000, thin=50, output=c("CI","S","TAU"))
#'@export
AM_mcmc_parameters <- function(  niter=5000,
                                 burnin=2500, ## niter / 2
                                 thin=1,
                                 verbose = 1,
                                 output=c("CI","K"),
                                 parallel=T,
                                 output_dir = NULL) {
  
  
  return (list(type="AM_MCMC_PARAMETERS", 
               niter=niter, burnin=burnin, thin=thin,
               verbose=verbose, output=output, parallel=parallel,
               output_dir=output_dir));
  
  
}


#################################################################################
##### Prior on Number of Components Configuration functions
#################################################################################


#' Generate a configuration object that contains a Point mass prior.
#' 
#' Generate a configuration object that assigns a Point mass prior to the number of mixture components.
#' This is the simplest option and it requires to specify a value \eqn{\widetilde{M}} 
#' such that \eqn{Pr(M=\widetilde{M}) =1}. 
#'
#'@param Mstar      Fixed value  \eqn{\widetilde{M}} for the number of components. 
#'@return list to be used as \code{mix_components_prior} argument for \code{\link{AM_mcmc_fit}}. 
#'
#'@keywords prior
#'@seealso \code{\link{AM_mcmc_fit}}
#'@export
#' 
#'@examples
#' 
#' AM_mix_components_prior_dirac (Mstar=3)
AM_mix_components_prior_dirac <- function(Mstar) {
  
  parameters = list(type = "AM_mix_components_prior_dirac", Mstar = Mstar);
  
  return (parameters);
};




#' Negative Binomial Prior.
#' 
#' This generate a configuration object for a Shifted Negative Binomial prior on the number of mixture components such as 
#'  \deqn{q_M(m)=Pr(M=m) =\frac{\Gamma(r+m-1)}{(m-1)!\Gamma(r)} p^{m-1}(1-p)^r, \quad m=1,2,3,\ldots}
#' The hyper-parameters \eqn{p\in (0,1)}  (probability of success) and \eqn{r>0} (size) can either be fixed using \code{r} and \code{p}
#' or assigned appropriate prior distributions. 
#' In the latter case, we assume \eqn{p \sim Beta(a_P,b_P)} and \eqn{r \sim  Gamma(a_R,b_R)}. In AntMAN we assume the following 
#' parametrization of the Gamma density: 
#' \deqn{p(x\mid a,b )= \frac{b^a x^{a-1}}{\Gamma(a)} \exp\{ -bx \}, \quad x>0  }
#' 
#' 
#' If no arguments are provided, the default is \eqn{r = 1 , a_P = 1, b_P = 1}.
#' 
#' Additionnaly, when init_R and init_P are no specified, there is default values : 
#' \eqn{init_R = 1} and \eqn{init_P = 0.5}
#'
#'@param a_R      The shape parameter \eqn{a}  of the \eqn{Gamma(a,b)} prior distribution for \eqn{r}.
#'@param b_R      The  rate parameter \eqn{b} of the \eqn{Gamma(a,b)} prior distribution for \eqn{r}.
#'@param init_R   The initial value of \eqn{r}, when specifying \code{a_R} and \code{b_R}.
#'@param a_P      The parameter\eqn{a}  of the \eqn{Beta(a,b)} prior distribution for \eqn{p}.
#'@param b_P      The parameter \eqn{b}  of the \eqn{Beta(a,b)} prior distribution for \eqn{p}.
#'@param init_P   The inivial  value of \eqn{p}, when specifying \code{a_P} and \code{b_P}.
#'@param R        It allows  to fix  \eqn{r} to a specific value.
#'@param P        It allows  to fix  \eqn{p} to a specific value.
#'
#'@return A configuration list to be used as \code{mix_components_prior} argument for \code{\link{AM_mcmc_fit}}. 
#'
#'
#'@keywords prior
#'@seealso \code{\link{AM_mcmc_fit}}
#'@export
#' 
#'@examples
#' 
#' AM_mix_components_prior_negbin (R=1, P=1)
#' AM_mix_components_prior_negbin ()

AM_mix_components_prior_negbin <- function(a_R = NULL, b_R = NULL, a_P = NULL, b_P = NULL, R = NULL, P = NULL, 
                                           init_R = NULL, init_P = NULL) {
  
  paradox_error_R = "Please note that you cannot specify a_R,b_R and R. R specifies a fixed value.";
  paradox_error_P = "Please note that you cannot specify a_P,b_P and P. P specifies a fixed value.";
  
  parameters = list(type = "AM_mix_components_prior_negbin");
  
       if (!is.null(a_R) & !is.null(b_R) & !is.null(init_R) &  is.null(R)) {parameters = append(parameters, list(a_R = a_R, b_R = b_R, init_R = init_R));}
  else if (!is.null(a_R) & !is.null(b_R) &  is.null(init_R) &  is.null(R)) {parameters = append(parameters, list(a_R = a_R, b_R = b_R,                ));}
  else if ( is.null(a_R) &  is.null(b_R) &  is.null(init_R) & !is.null(R)) {parameters = append(parameters, list(fixed_R = R));}
  else if ( is.null(a_R) &  is.null(b_R) &  is.null(init_R) &  is.null(R)) {parameters = append(parameters, list(fixed_R = 1));}
  else {stop ( paradox_error_R );}
  
       if (!is.null(a_P) & !is.null(b_P) & !is.null(init_P) &  is.null(P)) {parameters = append(parameters, list(a_P = a_P, b_P = b_P, init_P = init_P));}
  else if (!is.null(a_P) & !is.null(b_P) &  is.null(init_P) &  is.null(P)) {parameters = append(parameters, list(a_P = a_P, b_P = b_P,                ));}
  else if ( is.null(a_P) &  is.null(b_P) &  is.null(init_P) & !is.null(P)) {parameters = append(parameters, list(fixed_P = P));}
  else if ( is.null(a_P) &  is.null(b_P) &  is.null(init_P) &  is.null(P)) {parameters = append(parameters, list(a_P = 1, b_P = 1));}
  else {stop ( paradox_error_P );}
  
  return (parameters);
};




#' Generate a configuration object for a Poisson prior on the number of mixture components.
#' 
#' This function generates a configuration object for a Shifted Poisson prior on the number 
#' of mixture components such that  
#' \deqn{q_M(m)=     Pr (M=m)= \frac{e^{-\Lambda}\Lambda^{m-1} }{(m-1)!}    ,      \quad m=1,2,3,\ldots}
#' The hyper-parameter \eqn{\Lambda} can either be fixed using \code{Lambda} 
#' or assigned a \eqn{Gamma(a,b)} prior distribution with \code{a} and \code{b}.
#' In AntMAN we assume the following parametrization of the Gamma density: 
#' \deqn{p(x\mid a,b )= \frac{b^a x^{a-1}}{\Gamma(a)} \exp\{ -bx \}, \quad x>0  }
#' 
#' If no arguments are provided, the default is a prior distribution with \code{a = 1} and \code{b = 1}.
#'
#'@param a      The shape parameter \code{a} of  the \eqn{Gamma(a,b)} prior distribution.
#'@param b      The rate  parameter \code{b} of the \eqn{Gamma(a,b)} prior distribution.
#'@param init   The  initial value for \eqn{\Lambda}, when specifying \code{a} and \code{b}.
#'@param Lambda It allows to set the   hyper-parameter \eqn{\Lambda} to  fixed value.
#'
#'@return A configuration list to be used as \code{mix_components_prior} argument for \code{\link{AM_mcmc_fit}}. 
#'
#'@keywords prior
#'@seealso \code{\link{AM_mcmc_fit}}
#'@export
#' 
#'@examples
#' 
#' ## See \code{\link{AM_uninorm_mix_hyperparams}} example.
#' components_prior   = AM_mix_components_prior_pois (init=3,  a=1, b=1) 
#' 

AM_mix_components_prior_pois <- function(a=NULL, b=NULL, Lambda=NULL, init=NULL) {
  
  paradox_error = "Please note that you cannot specify a,b,init and Lambda. Lambda specifies a fixed value.";
  
  ### Default value ###
  parameters = list(type = "AM_mix_components_prior_pois",  a = 1, b = 1);
  
  if (!is.null(a) & !is.null(b)) {
    parameters = list(type = "AM_mix_components_prior_pois",  a = a, b = b);
    if (!is.null(init)) {
      parameters = list(type = "AM_mix_components_prior_pois",  a = a, b = b, init = init)
    };
    
    if (!is.null(Lambda)) {
      stop ( paradox_error );
    };
    
  } else if (!is.null(Lambda)) {
    parameters = list(type = "AM_mix_components_prior_pois", Lambda = Lambda);
    
    if (!is.null(a) | !is.null(b) | !is.null(init)) {
      stop ( paradox_error );
    }
    
  } else {
    if (!is.null(a) | !is.null(b) | !is.null(init)) {
      stop ( paradox_error );
    }
  };
  
  return (parameters);
};



#################################################################################
##### Configuration functions for the prior on the mixture weights
#################################################################################



#' Generate a configuration object to specify a prior on the hyper-parameter \eqn{\gamma} for the Dirichlet prior on the 
#' mixture weights. 
#' 
#' Generate a configuration object to specify a prior on the hyper-parameter \eqn{\gamma} for the Dirichlet prior on the 
#' mixture weights. 
#' We assume \eqn{\gamma \sim  Gamma(a,b)}. Alternatively we can fix \eqn{\gamma} to a specific value.
#' Default is \eqn{\gamma=1/N}, where N is the number of observations. 
#'In AntMAN we assume the following 
#' parametrization of the Gamma density: 
#' \deqn{p(x\mid a,b )= \frac{b^a x^{a-1}}{\Gamma(a)} \exp\{ -bx \}, \quad x>0  }
#' 
#'@param a      The shape parameter a of the Gamma prior
#'@param b      The rate parameter b of the Gamma prior
#'@param init   The init value for \eqn{\gamma}, when we assume \eqn{\gamma} random.
#'@param gamma  It allows to fix \eqn{\gamma}  to a specific value.
#'
#'@return A configuration list to be used as \code{mix_components_prior} argument for \code{\link{AM_mcmc_fit}}. 
#'
#'@examples 
#' AM_mix_weights_prior_gamma (a=1, b=1)
#' AM_mix_weights_prior_gamma (a=1, b=1, init=1)
#' AM_mix_weights_prior_gamma (gamma = 3)
#' AM_mix_weights_prior_gamma () 
#'@export
#'@keywords prior
AM_mix_weights_prior_gamma <- function(a = NULL, b = NULL, gamma = NULL, init = NULL) {
  
  paradox_error = "Please note that you cannot specify a,b,init and gamma. gamma specifies a fixed value.";
  
  parameters = list(type = "AM_mix_weights_prior_gamma");
  if (!is.null(a) & !is.null(b)) {
    parameters = list(type = "AM_mix_weights_prior_gamma",  a = a, b = b);
    if (!is.null(init)) {
      parameters = list(type = "AM_mix_weights_prior_gamma",  a = a, b = b, init = init)
    };
    
    if (!is.null(gamma)) {
      stop ( paradox_error );
    };
    
  } else if (!is.null(gamma)) {
    parameters = list(type = "AM_mix_weights_prior_gamma",  gamma = gamma);
    
    if (!is.null(a) | !is.null(b) | !is.null(init)) {
      stop ( paradox_error );
    }
    
  } else {
    if (!is.null(a) | !is.null(b) | !is.null(init)) {
      stop ( paradox_error );
    }
  };
  
  return (parameters);
};


#################################################################################
##### Mixture Kernel Hyperparameters.
#################################################################################


#'Univariate Poisson Mixture Hyperparameters.
#'
#' Generate a configuration object that specifies a univariate Poisson mixture kernel and allows to 
#' specify the hyperparameters of the  conjugate Gamma prior, i.e. the kernel is a \eqn{Poisson(\tau) }
#' and \eqn{\tau\sim Gamma(\alpha_0,\beta_0)}. 
#' In AntMAN we assume the following 
#' parametrization of the Gamma density: 
#' \deqn{p(x\mid a,b )= \frac{b^a x^{a-1}}{\Gamma(a)} \exp\{ -bx \}, \quad x>0  }
#'
#' Note, by default alpha0=1 and beta0=1.
#'
#'
#'
#'@param alpha0       The shape  hyperparameter \eqn{\alpha_0}.
#'@param beta0        The  rate hyperparameter \eqn{\beta_0}.
#'@return A list to be used as \code{mix_kernel_hyperparams} argument for \code{mcmc_fit}.
#'@examples 
#' AM_unipois_mix_hyperparams (alpha0=2, beta0=0.2)
#'@export
AM_unipois_mix_hyperparams <- function(alpha0, beta0) {
  return ( list ( type = "AM_unipois_mix_hyperparams", alpha0 = alpha0, beta0 = beta0 ) );
}

#' Univariate Normal Mixture Hyperparameters
#' 
#' Generate a configuration object that specifies univariate Normal mixture kernel and allows to set the hyperparameters of the Normal-InverseGamma conjugate prior. As such, the kernel is a Gaussian dsistribution 
#' with mean \eqn{\mu} and variance \eqn{\sigma^2}. The prior on \eqn{(\mu,\sigma^2)} the Normal-InverseGamma:
#' \deqn{\pi(\mu,\sigma^2\mid m_0,\kappa_0,\nu_0,\sigma^2_0) = \pi_{\mu}(\mu|\sigma^2,m_0,\kappa_0)\pi_{\sigma^2}(\sigma^2\mid \nu_0,\sigma^2_0)}
#'  \deqn{\pi_{\mu}(\mu|\sigma^2,m_0,\kappa_0)  =\frac{\sqrt{\kappa_0}}{\sqrt{2\pi\sigma^2}} 
#'  \exp^{-\frac{\kappa_0}{2\sigma^2}(\mu-m_0)^2 }, \qquad \mu\in\mathcal{R}}
#'  \deqn{\pi_{\sigma^2}(\sigma^2\mid \nu_0,\sigma^2_0)= {\frac {\sigma_0^{2^{\nu_0 }}}{\Gamma (\nu_0 )}}(1/\sigma^2)^{\nu_0 +1}\exp \left(-\frac{\sigma_0^2}{\sigma^2}\right), \qquad \sigma^2>0}
#' 
#' 
#' 
#' where \eqn{m_0} corresponds \code{m0}, 
#'       \eqn{\kappa_0} corresponds \code{k0}, 
#'       \eqn{\nu_0} corresponds \code{nu0}, 
#'       \eqn{\sigma^2_0} corresponds \code{sig02}. 
#' 
#'If hyperparameters are not specified, the default is \code{m0=0}, \code{k0=1}, \code{nu0=3}, \code{sig02=1}.
#'
#'@param m0      The \eqn{m_0} hyperparameter.
#'@param k0      The \eqn{\kappa_0} hyperparameter.
#'@param nu0     The \eqn{\nu_0} hyperparameter.
#'@param sig02   The \eqn{\sigma^2_0} hyperparameter.
#'@return A list to be used as \code{mix_kernel_hyperparams} argument for \code{mcmc_fit}. 
#'@examples 
#'      
#'      #### This example ...
#'      
#'      data(galaxy)
#'      y_uvn = galaxy
#'      mixture_uvn_params = AM_uninorm_mix_hyperparams  (m0=20.83146, k0=0.3333333,
#'                                                        nu0=4.222222, sig02=3.661027)
#'      
#'      mcmc_params        = AM_mcmc_parameters(niter=20000, burnin=5000, thin=10, verbose=1)
#'      components_prior   = AM_mix_components_prior_pois (init=3,  a=1, b=1) 
#'      weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)
#'      
#'      fit <- AM_mcmc_fit(
#'        y = y_uvn,
#'        mix_kernel_hyperparams = mixture_uvn_params,
#'        mix_components_prior =components_prior,
#'        mix_weight_prior = weights_prior,
#'        mcmc_parameters = mcmc_params)
#'      
#'      summary (fit)
#'      plot (fit)
#'      
#'@export
AM_uninorm_mix_hyperparams <- function(m0, k0, nu0, sig02) {
  return ( list ( type = "AM_uninorm_mix_hyperparams", m0 = m0 , k0 = k0 , nu0 = nu0 , sig02 = sig02 ) );
}

#' Univariate Binomial Mixture Hyperparameters.
#'  
#' Generate a configuration object that specifies the prior hyperparameters for a mixture of  Univariate Binomial kernels wth probability of success \eqn{\tau} and size \eqn{N}. 
#' The conjugate prior on \eqn{\tau} is a Beta distribution: 
#' \deqn{ \pi(\tau\mid \alpha,\beta)=\frac{\Gamma(\alpha+\beta)}{\Gamma(\alpha)\Gamma(\beta)} \tau^{\alpha-1}\left( 1-\tau \right)^{\beta-1} , \qquad 0\le\tau\le1}
#' \eqn{N} is fixed by the user and should always be specified. Here,  
#'  \eqn{\alpha} corresponds to \code{a0},\eqn{\beta} to \code{b0}.
#'  The default for the prior hyperparameters is \eqn{a0=1, b0=1}.
#' 
#' 
#'@param a0        The a0 hyperparameter.
#'@param b0        The b0 hyperparameter.
#'@param N         size of the Binomial distribution.
#'@return A list to be used as \code{mix_kernel_hyperparams} argument for \code{mcmc_fit}.
#'@examples 
#' AM_unibin_mix_hyperparams (a0=1,b0=1,N=100)
#'@export
AM_unibin_mix_hyperparams <- function(a0, b0, N) {
  ## TODO: Warning N > 1 an for small N almost not identifiable.			  
  return ( list ( type = "AM_unibin_mix_hyperparams",a0 = a0 , b0 = b0  , mb = N ) );
}


#' Multivariate Bernoulli Mixture Hyperparameters (Latent Class analysis)
#' 
#' Generate a configuration object that defines the prior hyperparameters for a mixture of multivariate Bernoulli.
#' If the dimension of the data is P, then the prior is a product of P independent Beta distributions, Beta(\eqn{a_{0i},a_{0i}}). Therefore,
#' the vectors of hyperparameters, a0 and b0,  are P-dimensional. Default is (a0= c(1,....,1),b0= c(1,....,1))
#'
#'@param a0        The a0 hyperparameters.
#'@param b0        The b0 hyperparameters.
#'@return A list to be used as \code{mix_kernel_hyperparams} argument for \code{mcmc_fit}.
#'@examples 
#' AM_multiber_mix_hyperparams (a0= c(1,1,1,1),b0= c(1,1,1,1))
#'@export
AM_multiber_mix_hyperparams <- function(a0, b0) {
  return ( list ( type = "AM_multiber_mix_hyperparams", a0 = a0 , b0 = b0  ) );
}




#' Multivariate Normal Mixture Hyperparameters.
#' 
#' 
#' This fnctions allows the user to specify the hyperparameters for the conjugate prior for a mixture of Multivariate Normals. We assume that the data are d-dimensional vectors \eqn{y_i}, where \eqn{y_i} are iid 
#' Normal randm variables with mean \eqn{\boldsymbol{\mu}} and covariance matrix \eqn{\boldsymbol{\Sigma}}.
#' The conjugate prior is 
#' \deqn{\pi(\boldsymbol \mu, \boldsymbol \Sigma\mid\boldsymbol m_0,\kappa_0,\nu_0,\boldsymbol \Lambda_0)= 
#' \pi_{\mu}(\boldsymbol \mu|\boldsymbol \Sigma,\boldsymbol m_0,\kappa_0)\pi_{\Sigma}(\boldsymbol \Sigma \mid \nu_0,\boldsymbol \Lambda_0)}
#'  \deqn{ \pi_{\mu}(\boldsymbol \mu|\boldsymbol \Sigma,\boldsymbol m_0,\kappa_0)  = 
#'  \frac{\sqrt{\kappa_0^d}}{\sqrt {(2\pi )^{d}|{\boldsymbol \Sigma }|}} \exp \left(-{\frac {\kappa_0}{2}}(\boldsymbol\mu -{\boldsymbol m_0 })^{\mathrm {T} }{\boldsymbol{\Sigma }}^{-1}(\boldsymbol\mu-{\boldsymbol m_0 })\right) \qquad \boldsymbol \mu\in\mathcal{R}^d}
#' \deqn{\pi_{\Sigma}(\boldsymbol \Sigma\mid \nu_0,\boldsymbol \Lambda_0)= {\frac {\left|{\boldsymbol \Lambda_0 }\right|^{\nu_0 /2}}{2^{\nu_0 d/2}\Gamma _{d}({\frac {\nu_0 }{2}})}}\left|\boldsymbol \Sigma \right|^{-(\nu_0 +d+1)/2}e^{-{\frac {1}{2}}\mathrm {tr} (\boldsymbol \Lambda_0 \boldsymbol \Sigma^{-1})}
#', \qquad \boldsymbol \Sigma^2>0}
#' with \code{mu0} corresponds to \eqn{\boldsymbol m_0}, \code{ka0} corresponds to  \eqn{\kappa_0}, 
#' \code{nu0} to \eqn{\nu_0}, \code{Lam0} to \eqn{\Lambda_0}.
#' 
#' Default is \code{(mu0=c(0,..,0)}, \code{ka0=1}, \code{nu0=Dim+2}, \code{Lam0=diag(Dim))} with \code{Dim} is the dimension of the data \code{y}.
#' We advise the user to set \eqn{\nu_0} equal to at least the dimension of the data, \code{Dim}, plus 2 
#'
#'@param mu0    The hyperparameter \eqn{\boldsymbol m_0}.
#'@param ka0    The hyperparameter \eqn{\kappa_0}.
#'@param nu0    The hyperparameter \eqn{\nu_0}.
#'@param Lam0   The hyperparameter \eqn{\Lambda_0}.
#'@return A list to be used as \code{mix_kernel_hyperparams} argument for \code{mcmc_fit}.
#'@examples 
#' AM_multinorm_mix_hyperparams ()
#'@export
AM_multinorm_mix_hyperparams <- function(mu0 = NULL, ka0 = NULL, nu0 = NULL, Lam0 = NULL) {
  return ( list ( type = "AM_multinorm_mix_hyperparams", mu0 = mu0 , ka0 = ka0 , nu0 = nu0 , Lam0 = Lam0 ) );
}

#################################################################################
##### Raffaele Functions
#################################################################################

#' Compute the logarithm of the absolute value of the  generalized Sriling number of second Kind (mi pare) See charambeloides, using a recursive formula Devo vedere la formula
#' 
#' There are no default values.
#'
#' @param n      The sample size
#' @param gamma    A positive real number \code{gamma} 
#'
#' @return A vector of length n, reporting the values \code{C(gamma,n,k)} for \code{k=1,...,n}
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' dd= AM_calcola_stirling_ricor_abs(11,10)
#' print(dd)
AM_calcola_stirling_ricor_abs <- function (n,gamma) {
	return(calcola_stirling_ricor_abs(n,gamma));
}

#' Compute the value V(n,k), needed to caclulate the eppf of a Finite Dirichlet process when the prior on the component-weigts of the mixture is a Dirichlet with parameter \code{gamma} (i.e. when unnormailized weights are distributed as Gamma(\eqn{\gamma},1) ) when the prior on the number of componet is Shifted Poisson of parameter \code{Lambda}. See Section 9.1.1 of the Paper Argiento de Iorio 2019.
#' 
#' There are no default values.
#'
#' @param n        The sample size
#' @param Lambda   The \code{Lambda} parameter of the Poisson
#' @param gamma    The \code{gamma} parameter of the Dirichlet 
#'
#' @return A vector of length n, reporting the values \code{V(n,k)} for \code{k=1,...,n}
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' n=1000
#' Lam=100
#' gam=0.5
#' vnk= AM_VnkPoisson(n,Lam,gam)
#' stir= AM_calcola_stirling_ricor_log(gam, n)
#' plot(exp(vnk+stir))
#' sum(exp(vvv+stir ))

AM_VnkPoisson <- function (n,Lambda,gamma) {
	return(VnkPoisson(n,Lambda,gamma));
}

#' This function compute the prior on the number of cluster, i.e. occupied component of the mixutre for a Finite Dirichlet process when the prior on the component-weigts of the mixture is a Dirichlet with parameter \code{gamma} (i.e. when unnormailized weights are distributed as Gamma(\eqn{\gamma},1) ) when the prior on the number of componet is Shifted Poisson of parameter \code{Lambda}. See Section 9.1.1 of Argiento de Iorio (2019) for more details.
#' 
#' There are no default values.
#'
#' @param n        The sample size
#' @param Lambda   The \code{Lambda} parameter of the Poisson
#' @param gamma    The \code{gamma} parameter of the Dirichlet 
#'
#' @return A vector of length n, reporting the values of the prior on the number of clusters induced by the prior on \code{M} and \code{w}, i.e. \code{p^*_k} for \code{k=1,...,n}. See Section 9.1.1 of Argiento de Iorio (2019) for more details.
#'
#' @keywords prior number of clusters
#'
#' @export
#' 
#' @examples
#' n <- 82
#' Lambda <- 10
#' gam_po <- 0.1550195
#' prior_K_po <-  AM_prior_K_Pois(n,gam_po,Lambda)
#' plot(1:n, prior_K_po, type = "n", bty = "l", xlab = "k", ylab = "P(K=k)",main="Prior on the number of clusters")

 AM_prior_K_Pois <- function (n,gamma,Lambda) {
  return(prior_K_Pois(n,gamma,Lambda));
 }


#' Compute the value V(n,k), needed to caclulate the eppf of a Finite Dirichlet process when the prior on the component-weigts of the mixture is a Dirichlet with parameter \code{gamma} (i.e. when unnormailized weights are distributed as Gamma(\eqn{\gamma},1) ) when the prior on the number of componet is Negative Binomial with parameter \code{r} and \code{p}with  mean is mu =1+ r*p/(1-p) TODO: CHECK THIS FORMULA!!!. See Section 9.1.1 of the Paper Argiento de Iorio 2019 for more details
#' 
#' There are no default values.
#'
#' @param n      The sample size
#' @param r      The dispersion parameter \code{r} of Negative Binomial
#' @param p      The probability of failure parameter \code{p} of Negative Binomial
#' @param gam    The \code{gamma} parameter of the Dirichlet 
#'
#' @return A vector of length n, reporting the values \code{V(n,k)} for \code{k=1,...,n}
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' n=1000
#' r=1000
#' p=0.5
#' gam=0.5
#' vnk= AM_VnkNegBin(n,r,p,gam);
#' stir= AM_calcola_stirling_ricor_log(gam, n)
#' plot(exp(vnk+stir+(1:n)*log(gam)))
#' sum(exp(vnk+stir))

AM_VnkNegBin <- function (n,r,p,gam) {
	return(VnkNegBin(n,r,p,gam));
}

#' This function compute the prior on the number of cluster, i.e. occupied component of the mixutre for a Finite Dirichlet process when the prior on the component-weigts of the mixture is a Dirichlet with parameter \code{gamma} (i.e. when unnormailized weights are distributed as Gamma(\eqn{\gamma},1) ) when the prior on the number of componet  is Negative Binomial with parameter \code{r>0} and \code{0<p<1}, with  mean is mu =1+ r*p/(1-p) TODO: CHECK THIS FORMULA!!!. See Section 9.1.1 of the Paper Argiento de Iorio 2019 for more details. 
#' 
#' There are no default values.
#'
#' @param n      The sample size
#' @param r      The dispersion parameter \code{r} of Negative Binomial
#' @param p      The probability of failure parameter \code{p} of Negative Binomial
#' @param gamma  The \code{gamma} parameter of the Dirichlet 
#'
#' @return A vector of length n, reporting the values \code{V(n,k)} for \code{k=1,...,n}
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' n <- 50
#' gamma <- 1
#' r <- 0.1
#' p <- 0.91
#' gam_nb <- 0.2381641
#' prior_K_nb <-  AM_prior_K_NegBin(n,gam_nb,r,p)
#' plot(1:n,prior_K_nb, type = "n", bty = "l", xlab = "k", ylab = "P(K=k)",main="Prior on the number of clusters")
#' lines(1:n,prior_K_de,type="h",lwd=2)

AM_prior_K_NegBin <- function (n,gam_nb, r, p){
	return(prior_K_NegBin(n,gam_nb, r, p));
}


#' Compute the value V(n,k), needed to caclulate the eppf of a Finite Dirichlet process when the prior on the component-weigts
#'  of the mixture is a Dirichlet with parameter \code{gamma} (i.e. when unnormailized weights are distributed as Gamma(\eqn{\gamma},1) )
#'   when the number of component are fixed to \code{M^*}, i.e. a Dirac prior assigning mass only to \code{M^*} is assumed. 
#'   See Section 9.1.1 of the Paper Argiento de Iorio 2019 for more details.
#' 
#' There are no default values.
#'
#' @param n      The sample size
#' @param Mstar  The number of component of the mixture 
#' @param gamma    The \code{gamma} parameter of the Dirichlet 
#'
#' @return A vector of length n, reporting the values \code{V(n,k)} for \code{k=1,...,n}
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' n=200
#' Mstar=100
#' gam=0.5
#' vvv=AM_VnkDelta(n,Mstar,gam);
#' stir= AM_calcola_stirling_ricor_log(gam, n)
#' stir
#' plot(exp(vvv+stir) )
#' sum(exp(vvv+stir ))

AM_VnkDelta <- function (n,Mstar,gam) {
	return(VnkDelta(n,Mstar,gam));
}

#' This function compute the prior on the number of cluster, i.e. occupied component of the mixutre for a Finite Dirichlet process 
#' when the prior on the component-weigts of the mixture is a Dirichlet with parameter \code{gamma} (i.e. when unnormailized weights 
#' are distributed as Gamma(\eqn{\gamma},1) ) when the number of component are fixed to \code{M^*}, i.e. a Dirac prior assigning mass
#'  only to \code{M^*} is assumed. See Section 9.1.1 of the Paper Argiento de Iorio 2019 for more details.#' There are no default values.
#'
#' @param n        The sample size
#' @param Mstar    The number of component of the mixture 
#' @param gamma    The \code{gamma} parameter of the Dirichlet 
#'
#' @return A vector of length n, reporting the values \code{V(n,k)} for \code{k=1,...,n}
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' n <- 82
#' gam_de <- 0.1743555
#' Mstar <- 12
#' prior_K_de <- AM_prior_K_Delta(n,gam_de, Mstar)
#' plot(1:n, prior_K_de, type = "n", bty = "l", xlab = "k", ylab = "P(K=k)",main="Prior on the number of clusters")

 AM_prior_K_Delta <- function (n,gam_de,Mstar){
 	return( prior_K_Delta(n,gam_de,Mstar));
  }


#' Once specified a fixed value of components \code{M^*} this function  adopt a  \emph{bisection method} to find the value of \code{gamma} 
#' such that the induced distribution on the number of clusers is centered around a user specifed value \eqn{K^*}, i.e. the function use
#'  a bisection method to solve Eq.~{eq:findgamma} of WE NEED TO CITE ANTMAN PAPER. The user can provide a lower \eqn{\gamma_{l}} and
#'  an upper \eqn{\gamma_{u}} bound for the possible values of $gamma$. The default values are \eqn{\gamma_l= 10^{-3}} and \eqn{\gamma_{u}=10}.
#'  A defaault value for the tolerance is \eqn{\epsilon=0.1}. Moreover, after a maximum number of iteration (default is 31), the function 
#'  stops warning that convergence has not bee reached.
#'
#' @param n             The sample size
#' @param Mstar         The number of component of the mixture 
#' @param Kstar         The mean number of cluster the user want to specify
#' @param gam_min=1e-4  The lower bound of the interval in which \code{gamma} should be lie
#' @param gam_max=10    The upper bound of the interval in which \code{gamma} should lie
#'
#'
#' @return A value of \code{gamma} such that \code{E(K)=K^*} 
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' n <- 82
#' Mstar <- 12
#' gam_de <- AM_find_gamma_Delta(n,Mstar,Kstar=6, gam_min=1e-4,gam_max=10, tolerance=0.1)
#' prior_K_de <-  AM_prior_K_Delta(n,gam_de,Mstar)
#' prior_K_de\%*\%1:n

AM_find_gamma_Delta <- function (n,Mstar,Kstar=6, gam_min=0.0001,gam_max=10, tolerance=0.1) {
	return(find_gamma_Delta(n,Mstar,Kstar, gam_min,gam_max, tolerance));
}


#' Once the prior on the numbuer of mixture $M$ is assumed to be a Shifted Posson of parameter \code{Lambda}, 
#' this function  adopt a \emph{bisection method} to find the value of \code{gamma} such that the induced distribution
#'  on the number of clusers is centered around a user specifed value \eqn{K^*}, i.e. the function use a bisection
#'   method to solve Eq.~{eq:findgamma} of WE NEED TO CITE ANTMAN PAPER. The user can provide a lower \eqn{\gamma_{l}} 
#'   and an upper \eqn{\gamma_{u}} bound for the possible values of $gamma$. The default values are \eqn{\gamma_l= 10^{-3}} and \eqn{\gamma_{u}=10}.
#'     A defaault value for the tolerance is \eqn{\epsilon=0.1}. Moreover, after a maximum number of iteration (default is 31),
#'      the function stops warning that convergence has not bee reached.
#'
#' @param n             The sample size
#' @param Lambda        The parameter of the Shifted Poisson for the number of components of the mixture
#' @param Kstar         The mean number of cluster the user want to specify
#' @param gam_min=1e-4  The lower bound of the interval in which \code{gamma} should be lie
#' @param gam_max=10    The upper bound of the interval in which \code{gamma} should lie
#'
#'
#' @return A value of \code{gamma} such that \code{E(K)=K^*} 
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' Lam  <- 11
#' gam_po <-  AM_find_gamma_Pois(n,Lam,Kstar=6, gam_min=0.0001,gam_max=10, tolerance=0.1)
#' prior_K_po <-  AM_prior_K_Pois(n,gam_po,Lam)
#' prior_K_po\%*\%1:n

 AM_find_gamma_Pois <- function (n,Lam,Kstar=6, gam_min=0.0001,gam_max=10, tolerance=0.1) {
 	return (find_gamma_Pois(n,Lam,Kstar, gam_min,gam_max, tolerance));
 }

#' Once the prior on the numbuer of mixture $M$ is assumed to be a Negative Binomial  Negative Binomial with parameter \code{r>0} and \code{0<p<1}, with  mean is 1+ r*p/(1-p), this function  adopt a \emph{bisection method} to find the value of \code{gamma} such that the induced distribution on the number of clusers is centered around a user specifed value \eqn{K^*}, i.e. the function use a bisection method to solve Eq.~{eq:findgamma} of WE NEED TO CITE ANTMAN PAPER. The user can provide a lower \eqn{\gamma_{l}} and an upper \eqn{\gamma_{u}} bound for the possible values of $gamma$. The default values are \eqn{\gamma_l= 10^{-3}} and \eqn{\gamma_{u}=10}.  A defaault value for the tolerance is \eqn{\epsilon=0.1}. Moreover, after a maximum number of iteration (default is 31), the function stops warning that convergence has not bee reached.
#'
#' @param n             The sample size
#' @param r      The dispersion parameter \code{r} of Negative Binomial
#' @param p      The probability of failure parameter \code{p} of Negative Binomial
#' @param Kstar         The mean number of cluster the user want to specify
#' @param gam_min=1e-4  The lower bound of the interval in which \code{gamma} should be lie
#' @param gam_max=10    The upper bound of the interval in which \code{gamma} should lie
#'
#'
#' @return A value of \code{gamma} such that \code{E(K)=K^*} 
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' r <- 1
#' p <- 0.8571
#' gam_nb= AM_find_gamma_NegBin(n,r,p,Kstar=6, gam_min=0.001,gam_max=10000, tolerance=0.1)
#' prior_K_nb= AM_prior_K_NegBin(n,gam_nb, r, p)
#' prior_K_nb\%*\%1:n

AM_find_gamma_NegBin <- function (n,r,p,Kstar=6, gam_min=0.001,gam_max=10000, tolerance=0.1){
	return (find_gamma_NegBin(n,r,p,Kstar, gam_min,gam_max, tolerance));
}


### This function takes as input the desired values of
### the marginal mean of mu = Emu
### the marginal variance of mu = Vmu
### the marginal mean of sig2 = Esig2 
### the marginal variance of sig2 = Vsig2
###
### The values returns a list with 
### hyperparameters of a Normal-Inverse-Gamma
### with the desired value of the marginal moments.
from_statitstics2hyperparam <- function(Emu,Vmu,Esig2,Vsig2){
	
	m0    = Emu
	nu0   = 2*(Esig2)^2/Vsig2+4
	sig02 = Esig2*(nu0-2)/nu0
	k0    = sig02/Vmu * nu0/(nu0-2)
	
	return(list(m0=m0,nu0=nu0,sig02=sig02,k0=k0))
}


#' this function compute the hyperparameters of an Normal-Inverse-Gamma 
#' distribution using an empirical Bayes approach.
#' @param  y       The data y
#' @param  scEmu a positive value scEmu (default =1) such that marginally E(mu)=(sample variance)*scEmu
#' @param  scEsig2 a positive value scEsig2 (default=3) such that marginally E(sig2)=(sample variance)*scEsig2
#' @param  CVsig2  The coefficient of variation of sig2 (default =3), CVsig2
#' @return Parameters of the normal inverse Gamma prior according to the empirical Bayes approach. 
#' @export
#' 

AM_emp_bayes_uninorm <- function(y,scEmu=1,scEsig2=3,CVsig2=3){
	n <- length(y)   ### sample size
	bary <- mean(y)  ### sample mean
	s2y <- var(y)    ### sample variance
	
	Emu <- bary
	Vmu <- s2y*scEmu
	Esig2 <- s2y/scEsig2
	Vsig2 <- CVsig2^2*Esig2^2
	
	return(from_statitstics2hyperparam(Emu,Vmu,Esig2,Vsig2))
	
}


