
#include "MixtureUnivariatePoisson.hpp"
#include "MixtureUnivariateNormal.hpp"


#include "MixtureMultiVariateNormal.hpp"
#include "MixtureUnivariateProbit.hpp"
#include "MixtureUnivariateBernoulli.hpp"
#include "MixtureMultivariateBernoulli.hpp"


#include "PriorPoisson.hpp"


static const  int AM_MIXTURE_UNIVARIATE            =  0;
static const  int AM_MIXTURE_MULTIVARIATE          =  1;

static const  int AM_MIXTURE_POISSON               =  2; // 001 [0/1]
static const  int AM_MIXTURE_BINOMIAL              =  4; // 010 [0/1]
static const  int AM_MIXTURE_NORMAL                =  6; // 011 [0/1]
static const  int AM_MIXTURE_PROBIT                =  8; // 100 [0/1]

static const  int AM_MIXTURE_UNIVARIATE_POISSON    =  AM_MIXTURE_UNIVARIATE | AM_MIXTURE_POISSON ;
static const  int AM_MIXTURE_UNIVARIATE_BINOMIAL   =  AM_MIXTURE_UNIVARIATE | AM_MIXTURE_BINOMIAL;
static const  int AM_MIXTURE_UNIVARIATE_NORMAL     =  AM_MIXTURE_UNIVARIATE | AM_MIXTURE_NORMAL  ;
static const  int AM_MIXTURE_UNIVARIATE_PROBIT     =  AM_MIXTURE_UNIVARIATE | AM_MIXTURE_PROBIT  ;

static const  int AM_MIXTURE_MULTIVARIATE_BINOMIAL = AM_MIXTURE_MULTIVARIATE | AM_MIXTURE_BINOMIAL;
static const  int AM_MIXTURE_MULTIVARIATE_NORMAL   = AM_MIXTURE_MULTIVARIATE | AM_MIXTURE_NORMAL;

constexpr bool is_univariate (int v) {return v & AM_MIXTURE_UNIVARIATE;};
constexpr bool is_multivariate (int v) {return v & AM_MIXTURE_MULTIVARIATE;};


static const int  AM_COMPONENTS_PRIOR_RANDOM_POISSON              = 1 ;
static const int  AM_COMPONENTS_PRIOR_FIXED_POISSON               = 2 ;
static const int  AM_COMPONENTS_PRIOR_FIXED_NEGATIVE_BINOMIAL     = 3 ;
static const int  AM_COMPONENTS_PRIOR_RANDOM_NEGATIVE_BINOMIAL    = 4 ;
static const int  AM_COMPONENTS_PRIOR_FIXED_DELTA_DIRAC           = 5 ;


static const int  AM_PRIOR_POISSON_GAMMA     = 1 ;
static const int  AM_PRIOR_GAMMA         	 = 2 ;
static const int  AM_PRIOR_FIXED         	 = 2 ;
static const int  AM_PRIOR_NEGATIVE_BINOMIAL = 3 ;

static const int  AM_GAMMA_WP       = 1 ;
static const int  AM_FIXED_WP       = 2 ;


// ************************************************ GIBBS ********************************************************************

// TODO facultavive initial clstering

Mixture* gen_mix (Rcpp::List          mixture_parameters) {
	int mixture_type = mixture_parameters["type"];
	Mixture*    mixture;
	switch(mixture_type) {
	case AM_MIXTURE_UNIVARIATE_POISSON :
		mixture = new Mixture_UnivariatePoisson (mixture_parameters["alpha0"], mixture_parameters["beta0"]);
		break;
	case AM_MIXTURE_UNIVARIATE_NORMAL :
		mixture = new Mixture_UnivariateNormal (mixture_parameters["m0"], mixture_parameters["k0"], mixture_parameters["nu0"], mixture_parameters["sig02"]);
		break;
	case AM_MIXTURE_UNIVARIATE_BINOMIAL :
		mixture = new Mixture_UnivariateBernoulli (mixture_parameters["a0"], mixture_parameters["b0"], mixture_parameters["mb"]) ;
		break;
	case AM_MIXTURE_UNIVARIATE_PROBIT :
		mixture = new Mixture_UnivariateProbit (mixture_parameters["X"], mixture_parameters["mu_beta"], mixture_parameters["Sig_beta"]) ;
		break;
	case AM_MIXTURE_MULTIVARIATE_NORMAL :
		mixture = new Mixture_MultivariateNormal (mixture_parameters["mu0"], mixture_parameters["ka0"], mixture_parameters["nu0"], mixture_parameters["Lam0"]);
		break;
	case AM_MIXTURE_MULTIVARIATE_BINOMIAL :
		mixture = new Mixture_MultivariateBernoulli (mixture_parameters["a0"], mixture_parameters["b0"], mixture_parameters["mb"]);
		break;
	default :
		mixture = NULL;
	}
	return mixture;
}

Prior* gen_prior (Rcpp::List mix_components_prior, Rcpp::List  mix_weight_prior ) {

	// Q is Fixed lambda or Poisson or Binomial
	// H is gamma

	Prior*                prior ;

	int prior_type   = mix_components_prior["type"];
	int weight_type  = mix_weight_prior["type"];

	assert(weight_type == AM_GAMMA_WP);
	double inith  = mix_weight_prior["init"];
	double ah  = mix_weight_prior["a"];
	double bh  = mix_weight_prior["b"];
	double lsd = mix_weight_prior["lsd"];



	switch (prior_type) {
	case AM_PRIOR_POISSON_GAMMA : {
		double aq    = mix_components_prior["a"];
		double bq    = mix_components_prior["b"];
		double initq = mix_components_prior["lambda"];
		prior = new PriorPoisson (
				inith,
				initq,
				ah, bh ,aq, bq, lsd);
		break;
	}
	case AM_PRIOR_GAMMA : {
		double initq = mix_components_prior["lambda"];
		prior = new PriorPoisson (inith,initq);
		break;
	}
	default :
		prior = NULL;
	}
	return prior;
}

// ************************************************ FUTURE ********************************************************************


// ************************************************ AM_mix_components_prior_ ********************************************************************

// [[Rcpp::export]]
Rcpp::List AM_mix_components_prior_pois(
		Rcpp::Nullable<double> a      = R_NilValue,
		Rcpp::Nullable<double> b      = R_NilValue,
		Rcpp::Nullable<double> Lambda = R_NilValue) { // TODO : can give only lambda or a and b

	static const std::string notice = "When using AM_mix_components_prior_pois, you need to specify either (lambda) or (a,b).";

	if (Lambda.isNull()) {

		VERBOSE_ASSERT(a.isNotNull() && b.isNotNull(), notice);
		 return Rcpp::List::create(
				 Rcpp::Named("type")     = AM_COMPONENTS_PRIOR_RANDOM_POISSON,
				 Rcpp::Named("a")        = a,
				 Rcpp::Named("b")        = b
				);
	} else {

		VERBOSE_ASSERT(a.isNull() && b.isNull(), notice);
		 return Rcpp::List::create(
				 Rcpp::Named("type")     = AM_COMPONENTS_PRIOR_FIXED_POISSON,
				 Rcpp::Named("lambda")   = Lambda
				);

	}


}


// [[Rcpp::export]]
Rcpp::List  AM_mix_components_prior_negbin(
		Rcpp::Nullable<double> a_R    = R_NilValue,
		Rcpp::Nullable<double> b_R    = R_NilValue,
		Rcpp::Nullable<double> a_P    = R_NilValue,
		Rcpp::Nullable<double> b_P    = R_NilValue,
		Rcpp::Nullable<double> LSDR_M = R_NilValue,
		Rcpp::Nullable<double> LSDP_M = R_NilValue,
		Rcpp::Nullable<double> R_M    = R_NilValue,
		Rcpp::Nullable<double> P_M    = R_NilValue){

	static const std::string notice = "When using AM_mix_components_prior_pois, you need to specify either (R_M,P_M) or (a,b,...).";

	if (R_M.isNotNull() and R_M.isNotNull()) {

		VERBOSE_ASSERT(a_R.isNull() && b_R.isNull(), notice);
		VERBOSE_ASSERT(a_P.isNull() && b_P.isNull(), notice);
		VERBOSE_ASSERT(LSDP_M.isNull() && LSDR_M.isNull(), notice);

		return Rcpp::List::create(
					 Rcpp::Named("type")        = AM_COMPONENTS_PRIOR_FIXED_NEGATIVE_BINOMIAL,
					 Rcpp::Named("R_M")   = R_M,
					 Rcpp::Named("P_M")   = P_M
					);
	}

	VERBOSE_ASSERT(R_M.isNull() && R_M.isNull(), notice);
	VERBOSE_ASSERT(a_R.isNotNull() && b_R.isNotNull(), notice);
	VERBOSE_ASSERT(a_P.isNotNull() && b_P.isNotNull(), notice);
	VERBOSE_ASSERT(LSDP_M.isNotNull() && LSDR_M.isNotNull(), notice);

	return Rcpp::List::create(
			 Rcpp::Named("type")  = AM_COMPONENTS_PRIOR_RANDOM_NEGATIVE_BINOMIAL,
			 Rcpp::Named("a_R")   = a_R,
			 Rcpp::Named("b_R")   = b_R,
			 Rcpp::Named("a_P")   = a_P,
			 Rcpp::Named("b_P")   = b_P,
			 Rcpp::Named("LSDR_M")   = LSDR_M,
			 Rcpp::Named("LSDP_M")   = LSDP_M
			);
} //TODO : lots of work to be done for different situations !! 4 combinations


// [[Rcpp::export]]
Rcpp::List AM_mix_components_prior_fixed(double Mstar){ // TODO: delta dirac !!
	 return Rcpp::List::create(
			 Rcpp::Named("type")        =  AM_COMPONENTS_PRIOR_FIXED_DELTA_DIRAC,
			 Rcpp::Named("Mstar")       = Mstar
			);
}

// ************************************************ AM_proposal_metropolis_hasting_params_ ********************************************************************

// TODO : lsd is step.

// ************************************************ AM_mix_weights_prior_ ********************************************************************


// [[Rcpp::export]]
Rcpp::List AM_mix_weights_prior_gamma(double a, double b, double lsd, double init) {
	 return Rcpp::List::create(
			 Rcpp::Named("type")   = AM_GAMMA_WP,
			 Rcpp::Named("a")   = a,
			 Rcpp::Named("b")   = b,
			 Rcpp::Named("lsd") = lsd,  // TODO: adaptive metropolis for 1.1 without LSD (no drugs in singapore!)
			 Rcpp::Named("init") = init // TODO: if not init, sample from the prior.
			);
}

// [[Rcpp::export]]
Rcpp::List AM_mix_weights_prior_fixed(double gamma) { // TODO: lsd not needed.
	 return Rcpp::List::create(
			 Rcpp::Named("type")   = AM_FIXED_WP,
			 Rcpp::Named("gamma")   = gamma
			);
}

// ************************************************ AM_mcmc_parameters ********************************************************************

// [[Rcpp::export]]
Rcpp::List AM_mcmc_parameters(
unsigned int niter=20000,
unsigned int burnin=10000,
unsigned int thin=10,
bool verbose=false,
std::string term_output = "all",
std::string file_output = "all",
std::string filename    = ""
) {
	 return Rcpp::List::create(
			 Rcpp::Named("niter")   = niter,
			 Rcpp::Named("burnin")   = burnin,
			 Rcpp::Named("thin")   = thin,
			 Rcpp::Named("verbose")   = verbose,
			 Rcpp::Named("term_output")   = term_output,
			 Rcpp::Named("file_output")   = file_output,
			 Rcpp::Named("filename")   = filename
				);
}


// ************************************************ mixture_uvp_hyperparams  ********************************************************************



		// [[Rcpp::export]]
		Rcpp::List AM_unipois_mix_hyperparams (double alpha0, double beta0) {

			return Rcpp::List::create(
					Rcpp::Named("alpha0")   = alpha0,
					Rcpp::Named("beta0")    = beta0,
					Rcpp::Named("type")     = AM_MIXTURE_UNIVARIATE_POISSON
					);

		}

		// [[Rcpp::export]]
		Rcpp::List AM_uninorm_mix_hyperparams (
				double m0, double k0, double nu0, double sig02
				) {
			return Rcpp::List::create(
					Rcpp::Named("m0")   = m0,
					Rcpp::Named("k0")   = k0,
					Rcpp::Named("nu0")   = nu0,
					Rcpp::Named("sig02")    = sig02,
					Rcpp::Named("type")     = AM_MIXTURE_UNIVARIATE_NORMAL
					);
		}

		// [[Rcpp::export]]
		Rcpp::List AM_unibin_mix_hyperparams (
				Rcpp::NumericVector a0,
				Rcpp::NumericVector b0,
				Rcpp::NumericVector mb
				) {
			return Rcpp::List::create(
					Rcpp::Named("a0")   = a0,
					Rcpp::Named("b0")   = b0,
					Rcpp::Named("mb")   = mb,
					Rcpp::Named("type")     = AM_MIXTURE_UNIVARIATE_BINOMIAL
					);
		}

		// [[Rcpp::export]]
		Rcpp::List AM_uniprobit_mix_hyperparams (
				const Rcpp::NumericMatrix & X,
				const Rcpp::NumericVector & Mu,
				const Rcpp::NumericMatrix & Sig
				) {
			return Rcpp::List::create(
					Rcpp::Named("X")   = X,
					Rcpp::Named("Mu")   = Mu,
					Rcpp::Named("Sig")   = Sig,
					Rcpp::Named("type")     = AM_MIXTURE_UNIVARIATE_PROBIT
					);
		}


		// [[Rcpp::export]]
		Rcpp::List AM_multibin_mix_hyperparams (
				const arma::vec  a0,
				const arma::vec  b0,
				const arma::vec  mb
				) {
			return Rcpp::List::create(
					Rcpp::Named("a0")   = a0,
					Rcpp::Named("b0")   = b0,
					Rcpp::Named("mb")   = mb,
					Rcpp::Named("type")     = AM_MIXTURE_MULTIVARIATE_BINOMIAL
					);
		}

		// [[Rcpp::export]]
		Rcpp::List AM_multinorm_mix_hyperparams (
				arma::vec    mu0,
				double       ka0,
				unsigned int nu0,
				arma::mat    Lam0) {
			return Rcpp::List::create(
					Rcpp::Named("mu0")   = mu0,
					Rcpp::Named("ka0")   = ka0,
					Rcpp::Named("nu0")   = nu0,
					Rcpp::Named("Lam0")   = Lam0,
					Rcpp::Named("type")     = AM_MIXTURE_MULTIVARIATE_NORMAL
					);
		}



// ************************************************ AM_mcmc_fit ********************************************************************

// [[Rcpp::export]]
Rcpp::List AM_mcmc_fit (
		Rcpp::RObject       y                      , /* Not optional */
		Rcpp::List          mix_kernel_hyperparams , /* Not optional */
		Rcpp::IntegerVector initial_clustering     = Rcpp::IntegerVector::create(), //  = Rcpp::IntegerVector::create() /* default will be 1for1 */
		int                 init_K                 = 0                            , //  = 0                             /* default will be 1for1 */
		Rcpp::List          mix_components_prior   = Rcpp::List::create()  , //  = Rcpp::List::create()          /* default will be Poisson (a,b)        */
		Rcpp::List          mix_weight_prior       = Rcpp::List::create()  , //  = Rcpp::List::create()          /* default will be (default  gamma = 1) */
		Rcpp::List          mcmc_parameters        = Rcpp::List::create()    //  = Rcpp::List::create()          /* (default niter=20000, ….) */
		) {


	// TODO: When you give me K, I need to do something,
	// TODO: When you give me CI, I up_ci next iteration
	// TODO: When you give me nothing, I init K = N, and up_ci second

	Prior*                prior    = gen_prior(mix_components_prior,mix_weight_prior);
	Mixture*              mixture  = gen_mix(mix_kernel_hyperparams);

	assert(mixture);
	assert(prior);

	int mixture_type = mix_kernel_hyperparams["type"];

	if(Rcpp::is<Rcpp::NumericVector>(y)){
	    assert (is_univariate(mixture_type)) ;
		GibbsResult res =  dynamic_cast<UnivariateMixture*>(mixture)->fit(Rcpp::as<arma::vec>(y) , Rcpp::as<cluster_indices_t>(initial_clustering), prior ,
				mcmc_parameters["niter"] ,mcmc_parameters["burnin"] ,mcmc_parameters["thin"] ,mcmc_parameters["verbose"] );

		return Rcpp::List::create(
					Rcpp::Named("U_post")      = res.U,
					Rcpp::Named("ci_post")     = res.ci,
					Rcpp::Named("S_post")      = res.S,
					Rcpp::Named("M_post")      = res.M,
					Rcpp::Named("K_post")      = res.K);

	} else if (Rcpp::is<Rcpp::NumericMatrix>(y)) {
		assert(is_multivariate(mixture_type));
		GibbsResult res =  dynamic_cast<MultivariateMixture*>(mixture)->fit(Rcpp::as<arma::mat>(y) , Rcpp::as<cluster_indices_t>(initial_clustering), prior ,
				mcmc_parameters["niter"] ,mcmc_parameters["burnin"] ,mcmc_parameters["thin"] ,mcmc_parameters["verbose"] );

		return Rcpp::List::create(
					Rcpp::Named("U_post")      = res.U,
					Rcpp::Named("ci_post")     = res.ci,
					Rcpp::Named("S_post")      = res.S,
					Rcpp::Named("M_post")      = res.M,
					Rcpp::Named("K_post")      = res.K);
	} else {
		VERBOSE_ERROR("Invalid arguments.");
	}


}
