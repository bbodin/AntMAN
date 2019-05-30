

##### Summary, Plot, and print functions for AntMan results

plot.AM_mcmc_fitness_result=function(x,...){
  if (!is.null(x$K)) plot(x$K,main="K Values")
  if (!is.null(x$CI)) plot(x$Y,col=x$CI[[length(x$CI)]]+1,main="Clusters")
}

summary.AM_mcmc_fitness_result=function(x,...){
	if (!is.null(x$K)) sprintf("Mean value of K is %f" , mean(x$K));
}


#' Generate a configuration object that contains parameters for a Poisson prior.
#' When there is no arguments, the default is a fixed Lambda = 1 / N, with N the input data size.
#'
#' @param a      The a parameter of the Poisson
#' @param b      The b parameter of the Poisson
#' @param init   The init value for Lambda of the Poisson
#' @param Lambda used to specify a fixed Lambda (instead of using a,b, and init).
#' @return A configuration list to be used as an argument for mcmc_fit. 
#' @examples 
#' AM_mix_components_prior_pois (a=1, b=1)
#' AM_mix_components_prior_pois (a=1, b=1, init=1)
#' AM_mix_components_prior_pois (Lambda = 3)
#' AM_mix_components_prior_pois () 
#' @export
AM_mix_components_prior_pois <- function(a = NULL, b = NULL, Lambda = NULL, init = NULL) {
  
  paradox_error = "Please note that you cannot specify a,b,init and Lambda. Lambda specifies a fixed value.";
  
  parameters = list(type = "AM_mix_components_prior_pois");
  
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
#' When there is no arguments, the default is *** TBD ***
#'
#' @param a_R      The a_R parameter of the Negative binomial
#' @param b_R      The b_R parameter of the Negative binomial
#' @param a_P      The a_R parameter of the Negative binomial
#' @param b_P      The b_R parameter of the Negative binomial
#' @param R_M      Used to specify a fixed R_M (instead of using a_R,b_R).
#' @param P_M      Used to specify a fixed P_M (instead of using a_P,b_P).
#' @return A configuration list to be used as an argument for mcmc_fit. 
#' @examples 
#' AM_mix_components_prior_negbin (M_R=1, M_P=1)
#' @export
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
#'
#' @param Mstar      Fixed value for M
#' @return A configuration list to be used as an argument for mcmc_fit. 
#' @examples 
#' AM_mix_components_prior_dirac (Mstar=3)
#' @export
AM_mix_components_prior_dirac <- function(Mstar) {
  
  parameters = list(type = "AM_mix_components_prior_dirac", Mstar = Mstar);
  
  return (parameters);
};



#' Generate a configuration object that contains parameters for a Gamma weight prior.
#'
#' @param a      The a parameter of the gamma
#' @param b      The b parameter of the gamma
#' @param init   The init value for Lambda of the gamma
#' @param gamma used to specify a fixed gamma (instead of using a,b, and init).
#' @return A configuration list to be used as an argument for mcmc_fit. 
#' @examples 
#' AM_mix_weights_prior_gamma (a=1, b=1)
#' AM_mix_weights_prior_gamma (a=1, b=1, init=1)
#' AM_mix_weights_prior_gamma (gamma = 3)
#' AM_mix_weights_prior_gamma () 
#' @export
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
#'
#' @param niter        Total number of iteration required.
#' @param burnin       Number of iteration to burn.
#' @param thin         Number of iteration to thin.
#' @param verbose      A value from 0 to 4, that specify the degres of verbosity (0:None,1:Warnings,2:Infos,4:Debug)
#' @param output       A list of output to return
#' @param file_output  A list of output to save in files
#' @param parallel     Some of the algorithms can be run in parallel using OpenMP. This parameter triggers the parallelism.
#' @return A configuration list to be used as an argument for mcmc_fit. 
#' @examples 
#' AM_mcmc_parameters (niter=1000, burnin=10000, thin=50)
#' AM_mcmc_parameters (niter=1000, burnin=10000, thin=50, output=c("CI","S","TAU"), file_output=c("all"))
#' @export
AM_mcmc_parameters <- function(  niter=20000,
                                 burnin=10000,
                                 thin=10,
                                 verbose = 1,
                                 output=c("CI","K"),
                                 parallel=0,
                                 file_output="M,K,Mna,Gamma") {
  
  
  return (list(type="AM_MCMC_PARAMETERS", 
            niter=niter, burnin=burnin, thin=thin,
            verbose=verbose, output=output, parallel=parallel,
            file_output=file_output));
  
  
}

#' Generate a configuration object that define univariate Poisson mixture hyperparameters.
#'
#' @param alpha0        The alpha0 hyperparameter.
#' @param beta0        The beta0 hyperparameter.
#' @return A configuration list to be used as an argument for mcmc_fit. 
#' @examples 
#' AM_unipois_mix_hyperparams (alpha0=2, beta0=0.2)
#' @export
AM_unipois_mix_hyperparams <- function(alpha0, beta0) {
  return ( list ( type = "AM_unipois_mix_hyperparams", alpha0 = alpha0, beta0 = beta0 ) );
}

#' Generate a configuration object that define univariate Normal mixture hyperparameters.
#'
#' @param m0      The m0 hyperparameter.
#' @param k0      The k0 hyperparameter.
#' @param nu0     The nu0 hyperparameter.
#' @param sig02   The sig02 hyperparameter.
#' @return A configuration list to be used as an argument for mcmc_fit. 
#' @examples 
#' AM_uninorm_mix_hyperparams (m0=0,k0=0.1,nu0=1,sig02=1.5)
#' @export
AM_uninorm_mix_hyperparams <- function(m0, k0, nu0, sig02) {
  return ( list ( type = "AM_uninorm_mix_hyperparams", m0 = m0 , k0 = k0 , nu0 = nu0 , sig02 = sig02 ) );
}

#' Generate a configuration object that define univariate binomial mixture hyperparameters.
#'
#' @param a0        The a0 hyperparameter.
#' @param b0        The b0 hyperparameter.
#' @param mb        The mb hyperparameter.
#' @return A configuration list to be used as an argument for mcmc_fit. 
#' @examples 
#' AM_unibin_mix_hyperparams (a0=1,b0=1,mb=100)
#' @export
AM_unibin_mix_hyperparams <- function(a0, b0, mb) {
  return ( list ( type = "AM_unibin_mix_hyperparams",a0 = a0 , b0 = b0  , mb = mb ) );
}


#' Generate a configuration object that define multivariate Bernoulli mixture hyperparameters.
#'
#' @param a0        The a0 hyperparameter.
#' @param b0        The b0 hyperparameter.
#' @return A configuration list to be used as an argument for mcmc_fit. 
#' @examples 
#' AM_multiber_mix_hyperparams (a0= c(1,1,1,1),b0= c(1,1,1,1))
#' @export
AM_multiber_mix_hyperparams <- function(a0, b0) {
  return ( list ( type = "AM_multiber_mix_hyperparams", a0 = a0 , b0 = b0  ) );
}

#' Generate a configuration object that define multivariate Normal mixture hyperparameters.
#'
#' @param mu0      The mu0 hyperparameter.
#' @param ka0      The ka0 hyperparameter.
#' @param nu0     The nu0 hyperparameter.
#' @param sig02   The sig02 hyperparameter.
#' @return A configuration list to be used as an argument for mcmc_fit. 
#' @examples 
#' AM_multinorm_mix_hyperparams (mu0=c(0,0),ka0=1,nu0=4,Lam0=diag(2))
#' @export
AM_multinorm_mix_hyperparams <- function(mu0, ka0, nu0, Lam0) {
  return ( list ( type = "AM_multinorm_mix_hyperparams", mu0 = mu0 , ka0 = ka0 , nu0 = nu0 , Lam0 = Lam0 ) );
}



#' Performs a Gibbs fit for the input data y, given a specific kernel (mix_kernel_hyperparams).
#'
#' @param y input data, can be a vector or a matrix.
#' @param mix_kernel_hyperparams is a configuration list, generated by *_mix_hyperparams functions.
#' @param initial_clustering is a vector CI of initial cluster assignement.
#' @param init_K is a prior on the number of cluster.
#' @param mix_components_prior is a configuration list generated with AM_mix_components_prior_* functions.
#' @param mix_weight_prior is a configuration list generated with AM_weight_prior_* functions.
#' @param mcmc_parameters is a configuration list generated with AM_mcmc_parameters. 
#' @return The output as specified by the mcmc parameters.
#' @examples
#' AM_mcmc_fit(AM_sample_unipois()$y, AM_unipois_mix_hyperparams (alpha0=2, beta0=0.2), mcmc_parameters = AM_mcmc_parameters(niter=2000, burnin=1000, thin=10))

 
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

