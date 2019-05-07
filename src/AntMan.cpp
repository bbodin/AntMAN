
#include "MixtureUnivariatePoisson.hpp"
#include "MixtureUnivariateNormal.hpp"


#include "MixtureMultiVariateNormal.hpp"
#include "MixtureUnivariateProbit.hpp"
#include "MixtureUnivariateBernoulli.hpp"
#include "MixtureMultivariateBernoulli.hpp"


/*
*/

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

static const int  AM_PRIOR_POISSON_GAMMA     = 1 ;
static const int  AM_PRIOR_RANDOM_POISSON    = 1 ;
static const int  AM_PRIOR_FIXED_POISSON     = 1 ;
static const int  AM_PRIOR_GAMMA         	 = 2 ;
static const int  AM_PRIOR_FIXED         	 = 2 ;
static const int  AM_PRIOR_NEGATIVE_BINOMIAL = 3 ;

static const int  AM_GAMMA_WP       = 1 ;
static const int  AM_FIXED_WP       = 2 ;


//TODO: optionnal parallel parameter to set openmp

// ********************************************* GENERAL *****************************************************************
// [[Rcpp::export]]
Rcpp::List AM_gibbs_parameters (
		int niter,
		int burnin,
		int thin,
		int verbose) {

	return Rcpp::List::create(
				Rcpp::Named("niter")   = niter,
				Rcpp::Named("burnin")  = burnin,
				Rcpp::Named("thin")    = thin,
				Rcpp::Named("verbose") = verbose);


}


// ************************************************ PRIOR PARAMETERS ********************************************************************

// [[Rcpp::export]]
Rcpp::List AM_poisson_gamma_prior_parameters (double inith , double initq, double ah, double bh,double aq, double bq,double lsd) {
	return Rcpp::List::create(
			Rcpp::Named("inith")= inith,
			Rcpp::Named("initq")= initq,
			Rcpp::Named("ah")   = ah,
			Rcpp::Named("bh")   = bh,
			Rcpp::Named("aq")   = aq,
			Rcpp::Named("bq")   = bq,
			Rcpp::Named("lsd")  = lsd,
			Rcpp::Named("type") = AM_PRIOR_POISSON_GAMMA
			);
}

// [[Rcpp::export]]
Rcpp::List AM_gamma_prior_parameters (double inith , double initq) {
	return Rcpp::List::create(
			Rcpp::Named("inith")= inith,
			Rcpp::Named("initq")= initq,
			Rcpp::Named("type") = AM_PRIOR_GAMMA
			);
}

// [[Rcpp::export]]
Rcpp::List AM_negative_binomial_prior_parameters (double initgamma, double R_M, double P_M , double LSDR_M, double LSDP_M,
	                                              double a_R, double b_R, double a_P, double b_P , double a, double b, double lsd) {

	return Rcpp::List::create(
			Rcpp::Named("initgamma")= initgamma,
			Rcpp::Named("R_M")    = R_M,
			Rcpp::Named("P_M")    = P_M,
			Rcpp::Named("LSDR_M") = LSDR_M,
			Rcpp::Named("LSDP_M") = LSDP_M,
			Rcpp::Named("a_R")    = a_R,
			Rcpp::Named("b_R")    = b_R,
			Rcpp::Named("a_P")    = a_P,
			Rcpp::Named("b_P")    = b_P,
			Rcpp::Named("a")      = a,
			Rcpp::Named("b")      = b,
			Rcpp::Named("lsd")    = lsd,
			Rcpp::Named("type")   = AM_PRIOR_NEGATIVE_BINOMIAL
			);
}


// ************************************************ MIXTURE PARAMETERS ********************************************************************


// [[Rcpp::export]]
Rcpp::List AM_univariate_poisson_mixture_parameters (double alpha0, double beta0) {

	return Rcpp::List::create(
			Rcpp::Named("alpha0")   = alpha0,
			Rcpp::Named("beta0")    = beta0,
			Rcpp::Named("type")     = AM_MIXTURE_UNIVARIATE_POISSON
			);

}

// [[Rcpp::export]]
Rcpp::List AM_univariate_normal_mixture_parameters (
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
Rcpp::List AM_univariate_binomial_mixture_parameters (
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
Rcpp::List AM_univariate_probit_mixture_parameters (
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
Rcpp::List AM_multivariate_binomial_mixture_parameters (
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
Rcpp::List AM_multivariate_normal_mixture_parameters (
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


// ************************************************ GIBBS ********************************************************************

// TODO facultavive initial clstering

Mixture* generate_mixture (Rcpp::List          mixture_parameters) {
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
Prior* generate_prior (Rcpp::List          prior_parameters) {

	int prior_type   = prior_parameters["type"];

	Prior*                prior ;
	switch (prior_type) {
	case AM_PRIOR_POISSON_GAMMA :
		prior = new PriorPoisson (prior_parameters["inith"], prior_parameters["initq"],
				prior_parameters["ah"], prior_parameters["bh"] , prior_parameters["aq"], prior_parameters["bq"], prior_parameters["lsd"]);
		break;
	case AM_PRIOR_GAMMA :
		prior = new PriorPoisson (prior_parameters["inith"], prior_parameters["initq"]);
		break;
	default :
		prior = NULL;
	}
	return prior;
}

// [[Rcpp::export]]
Rcpp::List AM_Univariate_Gibbs_Fit (
		Rcpp::NumericVector y,
		Rcpp::IntegerVector initial_clustering,
		Rcpp::List          prior_parameters,
		Rcpp::List          mixture_parameters,
		Rcpp::List          gibbs_parameters ) {

	Prior*                prior    = generate_prior(prior_parameters);
	UnivariateMixture*    mixture  = dynamic_cast<UnivariateMixture*>(generate_mixture(mixture_parameters));

	assert(mixture);
	assert(prior);

	GibbsResult res = mixture->fit(Rcpp::as<arma::vec>(y) , Rcpp::as<cluster_indices_t>(initial_clustering), prior ,
			   gibbs_parameters["niter"] ,gibbs_parameters["burnin"] ,gibbs_parameters["thin"] ,gibbs_parameters["verbose"] );

	return Rcpp::List::create(
				Rcpp::Named("U_post")      = res.U,
				Rcpp::Named("ci_post")     = res.ci,
				Rcpp::Named("S_post")      = res.S,
				Rcpp::Named("M_post")      = res.M,
				Rcpp::Named("K_post")      = res.K);

}

// [[Rcpp::export]]
Rcpp::List AM_Multivariate_Gibbs_Fit (
		Rcpp::NumericMatrix y,
		Rcpp::IntegerVector initial_clustering,
		Rcpp::List          prior_parameters,
		Rcpp::List          mixture_parameters,
		Rcpp::List          gibbs_parameters ) {

	Prior*                prior    = generate_prior(prior_parameters);
	MultivariateMixture*  mixture  = dynamic_cast<MultivariateMixture*>(generate_mixture(mixture_parameters));

	assert(mixture);
	assert(prior);

	GibbsResult res = mixture->fit(Rcpp::as<arma::mat>(y) , Rcpp::as<cluster_indices_t>(initial_clustering), prior ,
			   gibbs_parameters["niter"] ,gibbs_parameters["burnin"] ,gibbs_parameters["thin"] ,gibbs_parameters["verbose"] );

	return Rcpp::List::create(
				Rcpp::Named("U_post")      = res.U,
				Rcpp::Named("ci_post")     = res.ci,
				Rcpp::Named("S_post")      = res.S,
				Rcpp::Named("M_post")      = res.M,
				Rcpp::Named("K_post")      = res.K);

}


