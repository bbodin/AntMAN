
#################################################################################
##### Package Definition
#################################################################################

#' AntMAN: A package for fitting Finite Bayesian Mixture model with random number of component.
#'
#' AntMAN: Anthology of Mixture ANalysis tools 
#' AntMan is a R package to fit Finite Bayesian Mixture  model with random number of component. 
#' It is an alternative to Reversible Jump. 
#' Different kernels: Univariate Gaussian, Univariate Poisson, Univariate Binomial, Multivariate Gaussian, 
#' Multivariate Bernoulli (Latent Class Analysis). Conjugate priors.  
#'  We allow for different choices for the prior on the number of components: 
#' Shifted Poisson, Negative Binomial, and Point Masses (i.e. mixtures with fixed number of components).
#' 
#'@section Prior functions:
#' The Prior functions ...
#' 
#'@section Package Philosophy:
#' 
#' The main function of the AntMAN package is \code{AM_mcmc_fit} that performs a Gibbs sampling in order to estimate 
#' a mixture of a predifined type \code{mix_kernel_hyperparams} and that models a particular population given a sample \code{y}. 
#' Additionnaly a particular component prior \code{mix_components_prior} and a weight prior \code{mix_weight_prior} can be specified 
#' and \code{mcmc_parameters} will define the MCMC parameters for the Gibbs sampler (number of interation, burn-in, ...). 
#' A prior on the number of cluster (\code{init_K}) or a specific allocation (\code{init_clustering}) can also be specify. 
#' Otherwise, the default allocation we assign a different cluster for each element of the sample \code{y}.
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
#'@importFrom stats kmeans rbinom rnorm rpois runif 
#'@importFrom graphics plot
#'@description  
#'@docType package
#'@name AntMAN
NULL


#################################################################################
##### Datasets Definition
#################################################################################


#' Frequency of the word said in the Brown corpus.
#'     
#'@format  A list with 500 observations on the frequency of said in a text.
#'  
#'@source https://www.kaggle.com/nltkdata/brown-corpus
#'@references  Francis, W., and Kucera, H. (1982) Frequency Analysis of English Usage, Houghton Mifflin Company, Boston.
#'@examples
#'  data(said)
#'  
#'@keywords datasets
#'@docType data
"said"

#' galaxy
#'  
#' carcinoma is .... 
#'  
#'@format A data frame with X rows and Y variables.
#'  
#'@source URL?
#'  
#'@examples
#'  data(carcinoma)
#'  
#'@keywords datasets
#'@docType data
"galaxy"


#' carcinoma
#'  
#' carcinoma is .... 
#'  
#'@format A data frame with X rows and Y variables.
#'  
#'@source URL?
#'  
#'@examples
#'  data(carcinoma)
#'  
#'@keywords datasets
#'@docType data
"carcinoma"

#' brain
#'  
#' Teen Brain Images from the National Institutes of Health, U.S.
#'  
#' Picture of brain activity from a teenager consumming drugs. 
#'  
#'@format  A list of pixel taken from a YxY image in RGB format.
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
##### Summary, Plot, and print functions for AntMan results
#################################################################################

#' plot AM_mcmc_fitness_result  
#'  
#'  This fun...
#'  
#'@method plot AM_mcmc_fitness_result 
#'@export
plot.AM_mcmc_fitness_result=function(x,...){
  if (!is.null(x$K)) {
   hist(x$K,main="K values") ; 
   readline(prompt="Press [enter] to continue");
   }
  if (!is.null(x$M)) {hist(x$M,main="M values") ; readline(prompt="Press [enter] to continue");}
  if (!is.null(x$CI) && !is.null(x$Y)) {plot(x$Y,col=x$CI[[length(x$CI)]]+1,main="Clusters") ; readline(prompt="Press [enter] to continue");}
  if (!is.null(x$K)) {hist(x$K,main="Clusters") ; readline(prompt="Press [enter] to continue");}
  #### Histogram of  Gamma, ...
  #### Co clustering probability : How many time two are in the same custer. 
}

#'  summary AM_mcmc_fitness_result 
#'  
#'  This fun...
#'  
#'  
#'@method summary AM_mcmc_fitness_result 
#'@export
summary.AM_mcmc_fitness_result=function(object,...){
	print("Name\tMean\tStdDev") ;
	if (!is.null(object$K)) print(sprintf("%s\t%f\t%f","K" ,  mean(object$K) , sd(object$K))) ;
	if (!is.null(object$M)) print(sprintf("%s\t%f\t%f","M" ,  mean(object$M) , sd(object$M))) ;
	if (!is.null(object$Mna)) print(sprintf("%s\t%f\t%f","Mna" ,  mean(object$Mna) , sd(object$Mna))) ;
}


#################################################################################
##### Configuration functions
#################################################################################

#' Generate a configuration object that contains parameters for a Poisson prior.
#' 
#' When there is no arguments, the default is a  a=b=1
#'
#'@param a      The \code{a} parameter of the Poisson
#'@param b      The \code{b} parameter of the Poisson
#'@param init   The \code{init} value for \code{Lambda} of the Poisson
#'@param Lambda used to specify a fixed Lambda (instead of using \code{a},\code{b}, and \code{init}).
#'
#'@return output A configuration list to be used as an argument for mcmc_fit. 
#'
#'@keywords prior
#'
#'@export
#' 
#'@examples
#' AM_mix_components_prior_pois (a=1, b=1)
#' AM_mix_components_prior_pois (a=1, b=1, init=1)
#' AM_mix_components_prior_pois (Lambda = 3)
#' AM_mix_components_prior_pois () 
#' 

AM_mix_components_prior_pois <- function(a = NULL, b = NULL, Lambda = NULL, init = NULL) {
  
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



#' Generate a configuration object that contains parameters for a Negative Binomial prior.
#' When there is no arguments, the default is R = 1 , A_P =1 and B_P = 1
#'
#'@param a_R      The a_R parameter of the Negative binomial
#'@param b_R      The b_R parameter of the Negative binomial
#'@param a_P      The a_R parameter of the Negative binomial
#'@param b_P      The b_R parameter of the Negative binomial
#'@param R_M      Used to specify a fixed R_M (instead of using a_R,b_R).
#'@param P_M      Used to specify a fixed P_M (instead of using a_P,b_P).
#'@return A configuration list to be used as an argument for mcmc_fit. 
#'@examples 
#' AM_mix_components_prior_negbin (M_R=1, M_P=1)
#'@export
#'@keywords prior
AM_mix_components_prior_negbin <- function(a_R = NULL, b_R = NULL, a_P = NULL, b_P = NULL, M_R = NULL, M_P = NULL, init_R = NULL, init_P = NULL) {
  
  paradox_error_R = "Please note that you cannot specify a_R,b_R and R_M. R_M specifies a fixed value.";
  paradox_error_P = "Please note that you cannot specify a_P,b_P and P_M. P_M specifies a fixed value.";
  
  parameters = list(type = "AM_mix_components_prior_negbin");
   if (!is.null(a_R)) parameters = append(parameters, list(a_R = a_R));
   if (!is.null(b_R)) parameters = append(parameters, list(b_R = b_R));
   if (!is.null(init_R)) parameters = append(parameters, list(init_R = init_R));
   if (!is.null(M_R)) parameters = append(parameters, list(M_R = M_R));
   if (!is.null(a_P)) parameters = append(parameters, list(a_P = a_P));
   if (!is.null(b_P)) parameters = append(parameters, list(b_P = b_P));
   if (!is.null(init_P)) parameters = append(parameters, list(init_P = init_P));
   if (!is.null(M_P)) parameters = append(parameters, list(M_P = M_P));
    
  
  return (parameters);
};



#' Generate a configuration object that contains parameters for a Dirac prior.
#' No default value. 
#'
#'@param Mstar      Fixed value for M
#'@return A configuration list to be used as an argument for mcmc_fit. 
#'@examples 
#' AM_mix_components_prior_dirac (Mstar=3)
#'@export
#'@keywords prior
AM_mix_components_prior_dirac <- function(Mstar) {
  
  parameters = list(type = "AM_mix_components_prior_dirac", Mstar = Mstar);
  
  return (parameters);
};



#' Generate a configuration object that contains parameters for a Gamma weight prior.
#' Default Gamma = 1 / N
#'
#'@param a      The a parameter of the gamma
#'@param b      The b parameter of the gamma
#'@param init   The init value for Lambda of the gamma
#'@param gamma used to specify a fixed gamma (instead of using a,b, and init).
#'@return A configuration list to be used as an argument for mcmc_fit. 
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


#' Generate a configuration object that contains parameters for the MCMC.
#' burnin ===burnin  number of initial iteration to be dicarded, default niter/2  
#' niter= number of iteration after butnin  default 5000
#' thin, how oftern you save sample after burnn, i.e. one every thin, thin =1 save everything 
#' 
#' 
#' 
#' 
#' 
#'@param niter        Total number of iteration required.
#'@param burnin       Number of iteration to burn.
#'@param thin         Number of iteration to thin.
#'@param verbose      A value from 0 to 4, that specify the degres of verbosity (0:None,1:Warnings,2:Infos,4:Debug)
#'@param output       A list of output to return
#'@param file_output  A list of output to save in files
#'@param parallel     Some of the algorithms can be run in parallel using OpenMP. This parameter triggers the parallelism.
#'@return A configuration list to be used as an argument for mcmc_fit. 
#'@examples 
#' AM_mcmc_parameters (niter=1000, burnin=10000, thin=50)
#' AM_mcmc_parameters (niter=1000, burnin=10000, thin=50, output=c("CI","S","TAU"), file_output="")
#'@export
AM_mcmc_parameters <- function(  niter=20000,
                                 burnin=10000,
                                 thin=10,
                                 verbose = 1,
                                 output=c("CI","K"),
                                 parallel=0,
                                 file_output="") {
  
  
  return (list(type="AM_MCMC_PARAMETERS", 
            niter=niter, burnin=burnin, thin=thin,
            verbose=verbose, output=output, parallel=parallel,
            file_output=file_output));
  
  
}

#' Generate a configuration object that define univariate Poisson mixture hyperparameters.
#'alpha0=beta0=1
#'
#'
#'
#'@param alpha0        The alpha0 hyperparameter.
#'@param beta0        The beta0 hyperparameter.
#'@return A configuration list to be used as an argument for mcmc_fit. 
#'@examples 
#' AM_unipois_mix_hyperparams (alpha0=2, beta0=0.2)
#'@export
AM_unipois_mix_hyperparams <- function(alpha0, beta0) {
  return ( list ( type = "AM_unipois_mix_hyperparams", alpha0 = alpha0, beta0 = beta0 ) );
}

#' Generate a configuration object that define univariate Normal mixture hyperparameters.
#'default m0=0 k0=1 nu0= 3 sig02=1
#'
#'@param m0      The m0 hyperparameter.
#'@param k0      The k0 hyperparameter.
#'@param nu0     The nu0 hyperparameter.
#'@param sig02   The sig02 hyperparameter.
#'@return A configuration list to be used as an argument for mcmc_fit. 
#'@examples 
#'      
#'      #### This example ...
#'      
#'      data(galaxy)
#'      y_uvn = galaxy
#'      mixture_uvn_params = AM_uninorm_mix_hyperparams  (m0=20.83146, k0=0.3333333, nu0=4.222222, sig02=3.661027)
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

#' Generate a configuration object that define univariate binomial mixture hyperparameters.
#' a0=1 b0=1
#' 
#' 
#' 
#' 
#'@param a0        The a0 hyperparameter.
#'@param b0        The b0 hyperparameter.
#'@param mb        size of binomial .
#'@return A configuration list to be used as an argument for mcmc_fit. 
#'@examples 
#' AM_unibin_mix_hyperparams (a0=1,b0=1,mb=100)
#'@export
AM_unibin_mix_hyperparams <- function(a0, b0, mb) {
  return ( list ( type = "AM_unibin_mix_hyperparams",a0 = a0 , b0 = b0  , mb = mb ) );
}


#' Generate a configuration object that define multivariate Bernoulli mixture hyperparameters.
#'Default is (a0= c(1,....,1),b0= c(1,....,1))
#'
#'@param a0        The a0 hyperparameter.
#'@param b0        The b0 hyperparameter.
#'@return A configuration list to be used as an argument for mcmc_fit. 
#'@examples 
#' AM_multiber_mix_hyperparams (a0= c(1,1,1,1),b0= c(1,1,1,1))
#'@export
AM_multiber_mix_hyperparams <- function(a0, b0) {
  return ( list ( type = "AM_multiber_mix_hyperparams", a0 = a0 , b0 = b0  ) );
}

#' Generate a configuration object that define multivariate Normal mixture hyperparameters.
#' Default is (mu0=c(0,..,0),ka0=1,nu0=Dimension+2,Lam0=diag(Dim))
#'
#'
#'@param mu0      The mu0 hyperparameter.
#'@param ka0      The ka0 hyperparameter.
#'@param nu0     The nu0 hyperparameter.
#'@param sig02   The sig02 hyperparameter.
#'@return A configuration list to be used as an argument for mcmc_fit. 
#'@examples 
#' AM_multinorm_mix_hyperparams (mu0=c(0,0),ka0=1,nu0=4,Lam0=diag(2))
#'@export
AM_multinorm_mix_hyperparams <- function(mu0, ka0, nu0, Lam0) {
  return ( list ( type = "AM_multinorm_mix_hyperparams", mu0 = mu0 , ka0 = ka0 , nu0 = nu0 , Lam0 = Lam0 ) );
}

#################################################################################
##### AM_mcmc_fit function
#################################################################################


#' Performs a Gibbs sampling
#' 
#' The \code{AM_mcmc_fit} function performs a Gibbs sampling in order to estimate 
#' a mixture of a predifined type \code{mix_kernel_hyperparams} and that models a
#' particular population given a sample \code{y}. 
#' Additionnaly a particular component prior \code{mix_components_prior} and a weight 
#' prior \code{mix_weight_prior} can be specified and \code{mcmc_parameters} will define 
#' the MCMC parameters for the Gibbs sampler (number of interation, burn-in, ...). 
#' 
#' If no clustering is specified (either as \code{init_K} or \code{init_clustering}), 
#' then every observation is allocated a different clusters. 
#' If \code{init_K} is specified then we perform a Kmean. 
#' 
#' **Warning**: if you don't specify init_K or initial_cluster, the frist steps can be very long.
#'
#'@param y input data, can be a vector or a matrix.
#'@param mix_kernel_hyperparams is a configuration list, generated by *_mix_hyperparams functions.
#'@param initial_clustering is a vector CI of initial cluster assignement.
#'@param init_K is a prior on the number of cluster.
#'@param mix_components_prior is a configuration list generated with AM_mix_components_prior_* functions.
#'@param mix_weight_prior is a configuration list generated with AM_weight_prior_* functions.
#'@param mcmc_parameters is a configuration list generated with AM_mcmc_parameters. 
#'@return The return value is a \code{\link{AM_mcmc_fitness_result}} object. 
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
mix_components_prior = AM_mix_components_prior_pois() , 
mix_weight_prior = AM_mix_weights_prior_gamma(), 
mcmc_parameters = AM_mcmc_parameters() ) {
  
  if (is.null(init_K) & !is.null(initial_clustering)) {
  } else if (!is.null(init_K) & is.null(initial_clustering)) {
    initial_clustering <- kmeans(y, init_K)$cluster
  } else if (is.null(init_K) & is.null(initial_clustering)) {
    initial_clustering <- 0:(length(y)-1)
  } else {
    stop("Please provide either K_init or initial_clustering.")
  }
  
  structure(IAM_mcmc_fit(y = y, mix_kernel_hyperparams = mix_kernel_hyperparams, initial_clustering = initial_clustering, mix_components_prior = mix_components_prior, mix_weight_prior = mix_weight_prior, mcmc_parameters = mcmc_parameters)
            , class = "AM_mcmc_fitness_result") 
}

