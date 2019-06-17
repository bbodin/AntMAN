
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
#' @section Prior functions:
#' The Prior functions ...
#' The main function of the AntMAN package is \code{AM_mcmc_fit} that performs a Gibbs sampling in order to estimate a mixture of a predifined type \code{mix_kernel_hyperparams} and that models a particular population given a sample \code{y}. Additionnaly a particular component prior \code{mix_components_prior} and a weight prior \code{mix_weight_prior} can be specified and \code{mcmc_parameters} will define the MCMC parameters for the Gibbs sampler (number of interation, burn-in, ...). A prior on the number of cluster (\code{init_K}) or a specific allocation (\code{init_clustering}) can also be specify. Otherwise, the default allocation we assign a different cluster for each element of the sample \code{y}.
#' 
#' @section Package Philosophy:
#' For example, in order to identify clusters over a population of patients given a set of medical assumption:
#' 
#' mcmc = AM_mcmc_parameters(niter=20000) mix = AM_multiber_mix_hyperparams () fit = AM_mcmc_fit (mix, mcmc) summary (fit)
#' 
#' In this example \code{AM_multiber_mix_hyperparams} is one of the possible mixture to identify. AntMAN currently support five different mixtures :
#' 
#' AM_unipois_mix_hyperparams(alpha0, beta0) AM_uninorm_mix_hyperparams(m0, k0, nu0, sig02) AM_unibin_mix_hyperparams(a0, b0, mb) AM_multiber_mix_hyperparams(a0, b0) AM_multinorm_mix_hyperparams(mu0, ka0, nu0, Lam0)
#' 
#' Additionnaly, there is three prior_component available :
#' 
#' AM_mix_components_prior_pois AM_mix_components_prior_negbin AM_mix_components_prior_dirac
#' 
#' For example, in the context of image segmentation, where a maximal number of colour is require, a prior dirac can be used :
#' 
#' mcmc = AM_mcmc_parameters(niter=20000) mix = AM_multinorn_mix_hyperparams () prior_component = AM_mix_components_prior_dirac(10) # nothing more than 10 colours fit = AM_mcmc_fit (mix, prior_component, mcmc) summary (fit)
#' @importFrom Rcpp evalCpp
#' @importFrom stats kmeans rbinom rnorm rpois runif 
#' @importFrom graphics plot
#' @docType package
#' @name AntMAN
NULL


#################################################################################
##### Datasets
#################################################################################


#' said
#'  
#' Frequency of the word said in the Brown corpus.
#'  
#' This data set consider the frequency of appearance of the word said in the Borwn corpus.  
#'  
#'       For the record, here is how this data set has been generated : 
#'      
#'     	1. Download and unzip the data. Result is 500 files that contain series of sentences.
#'     
#'     >>> ls
#'     ca01  ca16  ca31  cb02  cb17  cc05  cd03  ce01  ce16  ce31  cf10  cf25  cf40  cg07  cg22  cg37  cg52  cg67  ch07  ch22  cj07  cj22  cj37  cj52  cj67  ck02  ck17  cl03  cl18  cn03  cn18  cp04  cp19  cr05
#'     ca02  ca17  ca32  cb03  cb18  cc06  cd04  ce02  ce17  ce32  cf11  cf26  cf41  cg08  cg23  cg38  cg53  cg68  ch08  ch23  cj08  cj23  cj38  cj53  cj68  ck03  ck18  cl04  cl19  cn04  cn19  cp05  cp20  cr06
#'     ca03  ca18  ca33  cb04  cb19  cc07  cd05  ce03  ce18  ce33  cf12  cf27  cf42  cg09  cg24  cg39  cg54  cg69  ch09  ch24  cj09  cj24  cj39  cj54  cj69  ck04  ck19  cl05  cl20  cn05  cn20  cp06  cp21  cr07
#'     ca04  ca19  ca34  cb05  cb20  cc08  cd06  ce04  ce19  ce34  cf13  cf28  cf43  cg10  cg25  cg40  cg55  cg70  ch10  ch25  cj10  cj25  cj40  cj55  cj70  ck05  ck20  cl06  cl21  cn06  cn21  cp07  cp22  cr08
#'     ca05  ca20  ca35  cb06  cb21  cc09  cd07  ce05  ce20  ce35  cf14  cf29  cf44  cg11  cg26  cg41  cg56  cg71  ch11  ch26  cj11  cj26  cj41  cj56  cj71  ck06  ck21  cl07  cl22  cn07  cn22  cp08  cp23  cr09
#'     ca06  ca21  ca36  cb07  cb22  cc10  cd08  ce06  ce21  ce36  cf15  cf30  cf45  cg12  cg27  cg42  cg57  cg72  ch12  ch27  cj12  cj27  cj42  cj57  cj72  ck07  ck22  cl08  cl23  cn08  cn23  cp09  cp24
#'     ca07  ca22  ca37  cb08  cb23  cc11  cd09  ce07  ce22  cf01  cf16  cf31  cf46  cg13  cg28  cg43  cg58  cg73  ch13  ch28  cj13  cj28  cj43  cj58  cj73  ck08  ck23  cl09  cl24  cn09  cn24  cp10  cp25
#'     ca08  ca23  ca38  cb09  cb24  cc12  cd10  ce08  ce23  cf02  cf17  cf32  cf47  cg14  cg29  cg44  cg59  cg74  ch14  ch29  cj14  cj29  cj44  cj59  cj74  ck09  ck24  cl10  cm01  cn10  cn25  cp11  cp26
#'     ca09  ca24  ca39  cb10  cb25  cc13  cd11  ce09  ce24  cf03  cf18  cf33  cf48  cg15  cg30  cg45  cg60  cg75  ch15  ch30  cj15  cj30  cj45  cj60  cj75  ck10  ck25  cl11  cm02  cn11  cn26  cp12  cp27
#'     ca10  ca25  ca40  cb11  cb26  cc14  cd12  ce10  ce25  cf04  cf19  cf34  cg01  cg16  cg31  cg46  cg61  ch01  ch16  cj01  cj16  cj31  cj46  cj61  cj76  ck11  ck26  cl12  cm03  cn12  cn27  cp13  cp28
#'     ca11  ca26  ca41  cb12  cb27  cc15  cd13  ce11  ce26  cf05  cf20  cf35  cg02  cg17  cg32  cg47  cg62  ch02  ch17  cj02  cj17  cj32  cj47  cj62  cj77  ck12  ck27  cl13  cm04  cn13  cn28  cp14  cp29
#'     ca12  ca27  ca42  cb13  cc01  cc16  cd14  ce12  ce27  cf06  cf21  cf36  cg03  cg18  cg33  cg48  cg63  ch03  ch18  cj03  cj18  cj33  cj48  cj63  cj78  ck13  ck28  cl14  cm05  cn14  cn29  cp15  cr01
#'     ca13  ca28  ca43  cb14  cc02  cc17  cd15  ce13  ce28  cf07  cf22  cf37  cg04  cg19  cg34  cg49  cg64  ch04  ch19  cj04  cj19  cj34  cj49  cj64  cj79  ck14  ck29  cl15  cm06  cn15  cp01  cp16  cr02
#'     ca14  ca29  ca44  cb15  cc03  cd01  cd16  ce14  ce29  cf08  cf23  cf38  cg05  cg20  cg35  cg50  cg65  ch05  ch20  cj05  cj20  cj35  cj50  cj65  cj80  ck15  cl01  cl16  cn01  cn16  cp02  cp17  cr03
#'     ca15  ca30  cb01  cb16  cc04  cd02  cd17  ce15  ce30  cf09  cf24  cf39  cg06  cg21  cg36  cg51  cg66  ch06  ch21  cj06  cj21  cj36  cj51  cj66  ck01  ck16  cl02  cl17  cn02  cn17  cp03  cp18  cr04
#'     
#'     
#'     
#'     	2. Analyse the shape of data, each sentences is a combinaison of words format like xxx/yyyy. 
#'     	To make sure of that, we replace every space by new line and makesure avery line has a '/' character.  
#'     
#'     >>> for f in c* ; do echo $f ;  cat $f | tr " " "\\n" | tr "\\t" "\\n" | grep -v "/" | sort -u ; done
#'     
#'     ca01
#'     
#'     ca02
#'     
#'     ca03
#'     
#'     ca04
#'     
#'     ....
#'     
#'     
#'     	3. for each file, we count the number of word.
#'     
#'     >>> cat ca01 | tr " " "\\n" | tr "\\t" "\\n" | grep -v  "^$" | sort -u| wc
#'         900     900    9706
#'     >>> cat ca01 | tr " " "\\n" | tr "\\t" "\\n" | grep -v  "^$" | wc
#'        2242    2242   20037
#'     
#'     
#'     	3bis. The average number of word is : 
#'     
#'     >>> for f in c* ; do cat $f | tr " " "\\n" | tr "\\t" "\\n" | grep -v  "^$" | awk ' { n++ } END { print n }'  ;  done | awk ' { sum += $1; n++ } END { print sum/n }' 
#'     2322.39
#'     
#'     	4. The word said, and it's different forms
#'     
#'     >>> for f in c* ; do cat $f | tr " " "\\n" | tr "\\t" "\\n" | grep -v  "^$" ;  done | grep "^[sS]aid/" | sort -u
#'     said/vbd
#'     Said/vbd
#'     said/vbn
#'     Said/vbn
#'     said/vbn-hl
#'     
#'     vbd is Verb, past tense
#'     vbn is Verb, past participle
#'     vbn-hl is Verb, past participle AND part of headline.
#'     
#'     	5. The frequency of the word "[Ss]aid/.*"
#'     
#'     >>> for f in c* ; do cat $f | tr " " "\\n" | tr "\\t" "\\n" | grep -v  "^$" | awk ' { n++ ; if ($1 == "said/vbd") m++ ;  if ($1 == "Said/vbd") m++ ;  if ($1 == "said/vbn") m++ ;  if ($1 == "Said/vbn") m++ ;  if ($1 == "said/vbn-hl") m++ ;  } END { print (2322.39 * m) / n; }'  ;  done >> said_frequency
#'     
#'     
#'  @format  A list with 500 observations on the frequency of said in a text.
#'  
#'  @source https://www.kaggle.com/nltkdata/brown-corpus
#'  @reference  Francis, W., and Kucera, H. (1982) Frequency Analysis of English Usage, Houghton Mifflin Company, Boston.
#'  @examples
#'  data(said)
#'  
#'  @keyword datasets
"said"

#' galaxy
#'  
#' carcinoma is .... 
#'  
#'  @format A data frame with X rows and Y variables.
#'  
#'  @source URL?
#'  
#'  @examples
#'  data(carcinoma)
#'  
#'  @keyword datasets
"galaxy"


#' carcinoma
#'  
#' carcinoma is .... 
#'  
#'  @format A data frame with X rows and Y variables.
#'  
#'  @source URL?
#'  
#'  @examples
#'  data(carcinoma)
#'  
#'  @keyword datasets
"carcinoma"

#' brain
#'  
#' Teen Brain Images from the National Institutes of Health, U.S.
#'  
#' Picture of brain activity from a teenager consumming drugs. 
#'  
#'  @format  A list of pixel taken from a YxY image in RGB format.
#'  
#'  @source https://www.flickr.com/photos/nida-nih/29741916012
#'  @references Crowley TJ, Dalwani MS, Mikulich-Gilbertson SK, Young SE, Sakai JT, Raymond KM, et al. (2015) Adolescents' Neural Processing of Risky Decisions: Effects of Sex and Behavioral Disinhibition. PLoS ONE 10(7): e0132322. doi:10.1371/journal.pone.0132322
#'  
#'  @examples
#'  data(brain)
#'  
#'  @keyword datasets
"brain"


#################################################################################
##### Summary, Plot, and print functions for AntMan results
#################################################################################

#' plot AM_mcmc_fitness_result  
#'  
#'  
#' @rdname AM_mcmc_fitness_result
#' @method plot AM_mcmc_fitness_result 
#' @export
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
#'  
#'  
#' @rdname AM_mcmc_fitness_result
#' @method summary AM_mcmc_fitness_result 
#' @export
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
#' @param a      The \code{a} parameter of the Poisson
#' @param b      The \code{b} parameter of the Poisson
#' @param init   The \code{init} value for \code{Lambda} of the Poisson
#' @param Lambda used to specify a fixed Lambda (instead of using \code{a},\code{b}, and \code{init}).
#'
#' @return output A configuration list to be used as an argument for mcmc_fit. 
#'
#' @keywords prior
#'
#' @export
#' 
#' @examples
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
#' No default value. 
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
#' Default Gamma = 1 / N
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
#' burnin ===burnin  number of initial iteration to be dicarded, default niter/2  
#' niter= number of iteration after butnin  default 5000
#' thin, how oftern you save sample after burnn, i.e. one every thin, thin =1 save everything 
#' 
#' 
#' 
#' 
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
#' AM_mcmc_parameters (niter=1000, burnin=10000, thin=50, output=c("CI","S","TAU"), file_output="")
#' @export
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
#'default m0=0 k0=1 nu0= 3 sig02=1
#'
#' @param m0      The m0 hyperparameter.
#' @param k0      The k0 hyperparameter.
#' @param nu0     The nu0 hyperparameter.
#' @param sig02   The sig02 hyperparameter.
#' @return A configuration list to be used as an argument for mcmc_fit. 
#' @examples 
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
#' @export
AM_uninorm_mix_hyperparams <- function(m0, k0, nu0, sig02) {
  return ( list ( type = "AM_uninorm_mix_hyperparams", m0 = m0 , k0 = k0 , nu0 = nu0 , sig02 = sig02 ) );
}

#' Generate a configuration object that define univariate binomial mixture hyperparameters.
#' a0=1 b0=1
#' 
#' 
#' 
#' 
#' @param a0        The a0 hyperparameter.
#' @param b0        The b0 hyperparameter.
#' @param mb        size of binomial .
#' @return A configuration list to be used as an argument for mcmc_fit. 
#' @examples 
#' AM_unibin_mix_hyperparams (a0=1,b0=1,mb=100)
#' @export
AM_unibin_mix_hyperparams <- function(a0, b0, mb) {
  return ( list ( type = "AM_unibin_mix_hyperparams",a0 = a0 , b0 = b0  , mb = mb ) );
}


#' Generate a configuration object that define multivariate Bernoulli mixture hyperparameters.
#'Default is (a0= c(1,....,1),b0= c(1,....,1))
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
#' Default is (mu0=c(0,..,0),ka0=1,nu0=Dimension+2,Lam0=diag(Dim))
#'
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

#################################################################################
##### Fit function
#################################################################################


#' Performs a Gibbs fit for the input data y, given a specific kernel (mix_kernel_hyperparams).
#' 
#' If no clustering is specified (either as K or init_cluster), then everybody in different clusters.
#' If K is specified then we perform a Kmean, 
#' 
#' Warning: if you don't specify init_K or initial_cluster, the frist steps can be very long.
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
#' @useDynLib AntMAN
#' @export
 
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

