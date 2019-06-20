
#################################################################################
##### Package Definition
#################################################################################

#' AntMAN: A package for fitting Finite Bayesian Mixture model with random number of component.
#'
#'@description AntMAN: Anthology of Mixture ANalysis tools 
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
#' The main function of the AntMAN package is \code{\link{AM_mcmc_fit}} that performs a Gibbs sampling in order to estimate 
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
#'@importFrom stats kmeans rbinom rnorm rpois runif sd
#'@importFrom graphics plot hist rasterImage
#'@importFrom sdols dlso
#'
#'
#'
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

#' S3 class AM_mcmc_output.
#' @description test
#' @exportClass AM_mcmc_output
#'@nameAM_mcmc_output
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
  if (!is.null(x$K)) {
   hist(x$K,main="K values") ; 
   readline(prompt="Press [enter] to continue");
   }
  if (!is.null(x$M)) {hist(x$M,main="M values") ; readline(prompt="Press [enter] to continue");}
  if (!is.null(x$K)) {hist(x$K,main="Clusters") ; readline(prompt="Press [enter] to continue");}
  #### Histogram of  Gamma, ...
  #### Co clustering probability : How many time two are in the same custer. 
  if (!is.null(x$CI) && !is.null(x$Y)) {plot(x$Y,col=x$CI[[length(x$CI)]]+1,main="Clusters") ; readline(prompt="Press [enter] to continue");}
}

#'  summary AM_mcmc_output 
#'  
#'  Print some useful informations about the mcmc results
#'  
#'@param object a AM_mcmc_output object
#'  
#'  
#'@method summary AM_mcmc_output 
#'@export
summary.AM_mcmc_output=function(object,...){
	print("Name\tMean\tStdDev") ;
	if (!is.null(object$K)) print(sprintf("%s\t%f\t%f","K" ,  mean(object$K) , sd(object$K))) ;
	if (!is.null(object$M)) print(sprintf("%s\t%f\t%f","M" ,  mean(object$M) , sd(object$M))) ;
	if (!is.null(object$Mna)) print(sprintf("%s\t%f\t%f","Mna" ,  mean(object$Mna) , sd(object$Mna))) ;
}


#################################################################################
##### AM_mcmc_fit function
#################################################################################


#' Performs a Gibbs sampling
#' 
#' The \code{AM_mcmc_fit} function performs a Gibbs sampling in order to estimate 
#' a mixture of a predifined type \code{mix_kernel_hyperparams} (generated with \code{AM_*_mix_hyperparams} functions) and that models a
#' particular population given a sample \code{y}. 
#' Additionnaly a particular component prior \code{mix_components_prior} (generated with  \code{AM_mix_components_prior_*} functions) and a weight 
#' prior \code{mix_weight_prior} (generated with  \code{AM_mix_weights_prior_*} functions) can be specified and \code{mcmc_parameters} will define 
#' the MCMC parameters forr the Gibbs sampler (generated with  \code{AM_mcmc_parameters} functions). 
#' 
#' If no clustering is specified (either as \code{init_K} or \code{init_clustering}), 
#' then every observation is allocated a different clusters. 
#' If \code{init_K} is specified then we perform a Kmean. 
#' 
#' **Warning**: if you don't specify init_K or initial_cluster, the fist steps can be very long because of default value.
#'
#'@param y input data, can be a vector or a matrix.
#'@param mix_kernel_hyperparams is a configuration list, generated by *_mix_hyperparams functions.
#'@param initial_clustering is a vector CI of initial cluster assignement.
#'@param init_K is a prior on the number of cluster.
#'@param mix_components_prior is a configuration list generated with AM_mix_components_prior_* functions.
#'@param mix_weight_prior is a configuration list generated with AM_weight_prior_* functions.
#'@param mcmc_parameters is a configuration list generated with AM_mcmc_parameters. 
#'@return The return value is a \code{AM_mcmc_output} object. 
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
            , class = "AM_mcmc_output") 
}

#################################################################################
##### AM_mcmc_parameters function
#################################################################################



#' MCMC Parameters
#' 
#' This generates an MCMC parameters list to be used as \code{mcmc_parameters} argument with \code{AM_mcmc_fit}. 
#' 
#' The \code{niter} argument specify the total number of iteration. 
#' \code{burnin} is the number of first iteration to drop.
#' \code{thin} specifies how oftern you save sample after burnin, i.e. one every thin, thin =1 save everything. 
#' 
#' 
#'@param niter        Total number of iteration required.
#'@param burnin       Number of iteration to burn.
#'@param thin         Thining rate.
#'@param verbose      A value from 0 to 4, that specify the degres of verbosity (0:None,1:Warnings,2:Infos,4:Debug)
#'@param output       A list of output to return
#'@param file_output  A list of output to save in files
#'@param parallel     Some of the algorithms can be run in parallel using OpenMP. This parameter triggers the parallelism.
#'@return list to be used as \code{mcmc_parameters} argument with \code{AM_mcmc_fit}. 
#'@examples 
#' AM_mcmc_parameters (niter=1000, burnin=10000, thin=50)
#' AM_mcmc_parameters (niter=1000, burnin=10000, thin=50, output=c("CI","S","TAU"), file_output="")
#'@export
AM_mcmc_parameters <- function(  niter=5000,
                                 burnin=2500, ## niter / 2
                                 thin=1,
                                 verbose = 1,
                                 output=c("CI","K"),
                                 parallel=0,
                                 file_output="") {
  
  
  return (list(type="AM_MCMC_PARAMETERS", 
               niter=niter, burnin=burnin, thin=thin,
               verbose=verbose, output=output, parallel=parallel,
               file_output=file_output));
  
  
}


#################################################################################
##### Prior Components Configuration functions
#################################################################################


#' Generate a configuration object that contains Point mass prior.
#' 
#' Generate a configuration object that contains Point mass prior.
#' This is the simplest option and it requires to specify a value \eqn{\widetilde{M}} 
#' such that \eqn{P(M=\widetilde{M}) =1}.
#'
#'@param Mstar      Fixed value for \eqn{\widetilde{M}} 
#'@return list to be used as \code{mix_components_prior} argument for \code{\link{AM_mcmc_fit}}. 
#'
#'@keywords prior
#'@seealso \code{\link{AM_mcmc_fit}}
#'@export
#' 
#'@examples
#' 
#' ## See \code{\link{???}} example.
#' AM_mix_components_prior_dirac (Mstar=3)
AM_mix_components_prior_dirac <- function(Mstar) {
  
  parameters = list(type = "AM_mix_components_prior_dirac", Mstar = Mstar);
  
  return (parameters);
};




#' Negative Binomial Prior.
#' 
#' This generate a configuration object for a Negative Binomial prior such as 
#'  \deqn{q_M(m)=\frac{\Gamma(r+m-1)}{(m-1)!\Gamma(r)} p^{m-1}(1-p)^r, \quad m=1,2,3,\ldots}
#' The hyper-parameters \eqn{p\in (0,1)} and \eqn{r>0} can either be fixed using \code{r} and \code{p}
#' or assigned prior distributions. 
#' In the latter case, we assume \eqn{p \sim Beta(a_P,b_P)} and \eqn{r \sim  Gamma(a_R,b_R)}
#' 
#' If no arguments are provided, the default is \eqn{r = 1 , a_P = 1 and b_P = 1}.
#'
#'@param a_R      The \eqn{a} parameter of \eqn{Gamma(a,b)} prior distribution for \eqn{r}.
#'@param b_R      The \eqn{b} parameter of \eqn{Gamma(a,b)} prior distribution for \eqn{r}.
#'@param init_R   The first value of \eqn{r} when using \code{a_R} and \code{b_R}.
#'@param a_P      The \eqn{a} parameter of \eqn{Beta(a,b)} prior distribution for \eqn{p}.
#'@param b_P      The \eqn{a} parameter of \eqn{Beta(a,b)} prior distribution for \eqn{p}.
#'@param init_P   The first value of \eqn{p} when using \code{a_P} and \code{b_P}.
#'@param R      Used to specify a fixed value \eqn{r}.
#'@param P      Used to specify a fixed value \eqn{p}.
#'@return A configuration list to be used as \code{mix_components_prior} argument for \code{\link{AM_mcmc_fit}}. 
#'
#'@keywords prior
#'@seealso \code{\link{AM_mcmc_fit}}
#'@export
#' 
#'@examples
#' 
#' ## See \code{\link{???}} example.
#' AM_mix_components_prior_negbin (R=1, P=1)

AM_mix_components_prior_negbin <- function(a_R = NULL, b_R = NULL, a_P = NULL, b_P = NULL, R = NULL, P = NULL, 
                                           init_R = NULL, init_P = NULL) {
  
  paradox_error_R = "Please note that you cannot specify a_R,b_R and R_M. R_M specifies a fixed value.";
  paradox_error_P = "Please note that you cannot specify a_P,b_P and P_M. P_M specifies a fixed value.";
  
  parameters = list(type = "AM_mix_components_prior_negbin");
  if (!is.null(a_R)) parameters = append(parameters, list(a_R = a_R));
  if (!is.null(b_R)) parameters = append(parameters, list(b_R = b_R));
  if (!is.null(init_R)) parameters = append(parameters, list(init_R = init_R));
  if (!is.null(R)) parameters = append(parameters, list(M_R = R));
  if (!is.null(a_P)) parameters = append(parameters, list(a_P = a_P));
  if (!is.null(b_P)) parameters = append(parameters, list(b_P = b_P));
  if (!is.null(init_P)) parameters = append(parameters, list(init_P = init_P));
  if (!is.null(P)) parameters = append(parameters, list(M_P = P));
  
  
  return (parameters);
};




#' Generate a configuration object for a Poisson prior.
#' 
#' This generate a configuration object for a Poisson prior such as 
#' \deqn{q_M(m)=     \frac{e^{-\Lambda}\Lambda^{m-1} }{(m-1)!}    ,      \quad m=1,2,3,\ldots}
#' The hyper-parameter \eqn{\Lambda} can either be fixed using \code{Lambda} 
#' or assigned a \eqn{Gamma(a,b)} prior distribution with \code{a} and \code{b}.
#' 
#' If no arguments are provided, the default is prior distribution with \code{a = 1} and \code{b = 1}.
#'
#'@param a      The \code{a} parameter of \eqn{Gamma(a,b)} prior distribution.
#'@param b      The \code{b} parameter of \eqn{Gamma(a,b)} prior distribution.
#'@param init   The \code{init} value for \eqn{\Lambda} when using \code{a} and \code{b}.
#'@param Lambda used to specify a fixed  hyper-parameter \eqn{\Lambda}.
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
##### Weights Prior Configuration functions
#################################################################################



#' Generate a configuration object that contains parameters for a Gamma weight prior.
#' Default Gamma = 1 / N
#'
#'@param a      The a parameter of the gamma
#'@param b      The b parameter of the gamma
#'@param init   The init value for Lambda of the gamma
#'@param gamma used to specify a fixed gamma (instead of using a,b, and init).
#'@return A configuration list to be used as \code{mix_weight_prior} argument for \code{mcmc_fit}. 
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
##### Mixture Hyperparameters.
#################################################################################


#' Univariate Poisson Mixture Hyperparameters.
#' 
#'  Generate a configuration object that define univariate Poisson mixture hyperparameters suc has :
#'  \deqn{\pi(\tau \mid \alpha,\beta)=\frac{\beta^\alpha}{\Gamma(\alpha)}\tau^{\alpha-1}\mathrm{e}^{-\beta\tau}, \qquad \tau >0}
#'
#'  By default alpha0=1 and beta0=1.
#'
#'
#'
#'@param alpha0        The alpha0 hyperparameter.
#'@param beta0        The beta0 hyperparameter.
#'@return A list to be used as \code{mix_kernel_hyperparams} argument for \code{mcmc_fit}.
#'@examples 
#' AM_unipois_mix_hyperparams (alpha0=2, beta0=0.2)
#'@export
AM_unipois_mix_hyperparams <- function(alpha0, beta0) {
  return ( list ( type = "AM_unipois_mix_hyperparams", alpha0 = alpha0, beta0 = beta0 ) );
}

#' Univariate Normal Mixture Hyperparameters
#' 
#' 
#' Generate a configuration object that define univariate Normal mixture hyperparameters as
#' \deqn{\pi(\mu,\sigma^2\mid m_0,\kappa_0,\nu_0,\sigma^2_0) = \pi_{\mu}(\mu|\sigma^2,m_0,\kappa_0)\pi_{\sigma^2}(\sigma^2\mid \nu_0,\sigma^2_0)}
#' Where \eqn{m_0} is \code{m0}, 
#'       \eqn{\kappa_0} is \code{k0}, 
#'       \eqn{\nu_0} is \code{nu0}, 
#'       \eqn{\sigma^2_0} is \code{sig02}. 
#' 
#'If hyperparameters are not, default is \code{m0=0}, \code{k0=1}, \code{nu0=3}, \code{sig02=1}.
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
#' Generate a configuration object that define Univariate Binomial Mixture Hyperparameters such as :
#' \deqn{ \pi(\theta\mid \alpha,\beta)=\frac{\Gamma(\alpha+\beta)}{\Gamma(\alpha)\Gamma(\beta)} \tau^{\alpha-1}\left( 1-\tau \right)^{\beta-1} , \qquad 0\le\tau\le1}
#'with \eqn{\alpha} is \code{a0},\eqn{\beta} is \code{b0}, and N is ... oups I need to add more infos in those sections.
#' a0=1 b0=1
#' 
#' 
#'@param a0        The a0 hyperparameter.
#'@param b0        The b0 hyperparameter.
#'@param mb        size of binomial .
#'@return A list to be used as \code{mix_kernel_hyperparams} argument for \code{mcmc_fit}.
#'@examples 
#' AM_unibin_mix_hyperparams (a0=1,b0=1,mb=100)
#'@export
AM_unibin_mix_hyperparams <- function(a0, b0, mb) {
  return ( list ( type = "AM_unibin_mix_hyperparams",a0 = a0 , b0 = b0  , mb = mb ) );
}


#' Multivariate Bernoulli Mixture Hyperparameters (Latent Class analysis)
#' 
#' Generate a configuration object that define multivariate Bernoulli mixture hyperparameters.
#'Default is (a0= c(1,....,1),b0= c(1,....,1))
#'
#'@param a0        The a0 hyperparameter.
#'@param b0        The b0 hyperparameter.
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
#' Generate Multivariate Normal Mixture hyperparameters such as 
#' \deqn{\pi(\boldsymbol \mu, \boldsymbol \Sigma\mid\boldsymbol m_0,\kappa_0,\nu_0,\boldsymbol \Lambda_0)= \pi_{\mu}(\boldsymbol \mu|\boldsymbol \Sigma,\boldsymbol m_0,\kappa_0)\pi_{\Sigma}(\boldsymbol \Sigma \mid \nu_0,\boldsymbol \Lambda_0)}
#' with \code{mu0} is \eqn{\boldsymbol m_0}, \code{ka0} is \eqn{\kappa_0}, 
#' \code{nu0} is \eqn{\nu_0}, \code{Lam0} is \eqn{\Lambda_0}.
#' 
#' Default is \code{(mu0=c(0,..,0),ka0=1,nu0=Dim+2,Lam0=diag(Dim))} with \code{Dim} the number of variable per observation.
#'
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

