
#include "MixtureUnivariatePoisson.hpp"
#include "MixtureUnivariateNormal.hpp"


#include "MixtureMultiVariateNormal.hpp"
#include "MixtureUnivariateProbit.hpp"
#include "MixtureUnivariateBernoulli.hpp"
#include "MixtureMultivariateBernoulli.hpp"


#include "PriorPoisson.hpp"

#include "utils.hpp"

bool is_univariate (Rcpp::List mix_kernel_hyperparams) {
	VERBOSE_ASSERT (mix_kernel_hyperparams.containsElementNamed("type"), "mix_kernel_hyperparams does not contain a type field.");
	std::string mixture_type = mix_kernel_hyperparams["type"];
	return mixture_type.find("uni") != std::string::npos;
};
bool is_multivariate (Rcpp::List mix_kernel_hyperparams) {
	VERBOSE_ASSERT (mix_kernel_hyperparams.containsElementNamed("type"), "mix_kernel_hyperparams does not contain a type field.");
	std::string  mixture_type = mix_kernel_hyperparams["type"];
	return mixture_type.find("multi") != std::string::npos;
};

static const int  AM_OUTPUT_CI  = 1 << 0;
static const int  AM_OUTPUT_TAU = 1 << 1;
static const int  AM_OUTPUT_S   = 1 << 2;
static const int  AM_OUTPUT_M   = 1 << 3;
static const int  AM_OUTPUT_K   = 1 << 4;
static const int  AM_OUTPUT_Mna = 1 << 5;
static const int  AM_OUTPUT_H   = 1 << 6;
static const int  AM_OUTPUT_Q   = 1 << 7;

static const int  AM_OUTPUT_DEFAULT   = AM_OUTPUT_CI & AM_OUTPUT_TAU & AM_OUTPUT_S & AM_OUTPUT_M & AM_OUTPUT_K;



// ************************************************ AM_mcmc_fit ********************************************************************

		Mixture* gen_mix (Rcpp::List          mixture_parameters) {

			VERBOSE_ASSERT (mixture_parameters.containsElementNamed("type"), "In gen_mix mixture_parameters does not contain a type field.");

			if (mixture_parameters.containsElementNamed("alpha0") and mixture_parameters.containsElementNamed("beta0")) {
				return new Mixture_UnivariatePoisson (mixture_parameters["alpha0"], mixture_parameters["beta0"]);
			}

			if (mixture_parameters.containsElementNamed("m0") and mixture_parameters.containsElementNamed("k0")
					 and mixture_parameters.containsElementNamed("nu0")
					 and mixture_parameters.containsElementNamed("sig02")) {
				return new Mixture_UnivariateNormal (mixture_parameters["m0"], mixture_parameters["k0"], mixture_parameters["nu0"], mixture_parameters["sig02"]);
			}

			if (mixture_parameters.containsElementNamed("a0") and mixture_parameters.containsElementNamed("b0")
							 and mixture_parameters.containsElementNamed("mb")) {
				return  new Mixture_UnivariateBernoulli (mixture_parameters["a0"], mixture_parameters["b0"], mixture_parameters["mb"]) ;
			}

			if (mixture_parameters.containsElementNamed("X") and mixture_parameters.containsElementNamed("mu_beta")
							 and mixture_parameters.containsElementNamed("Sig_beta")) {
				return new Mixture_UnivariateProbit (mixture_parameters["X"], mixture_parameters["mu_beta"], mixture_parameters["Sig_beta"]) ;
			}

			if (mixture_parameters.containsElementNamed("mu0") and mixture_parameters.containsElementNamed("ka0")
							 and mixture_parameters.containsElementNamed("nu0")
							 and mixture_parameters.containsElementNamed("Lam0")) {
				return new Mixture_MultivariateNormal (mixture_parameters["mu0"], mixture_parameters["ka0"], mixture_parameters["nu0"], mixture_parameters["Lam0"]);
			}

			if (mixture_parameters.containsElementNamed("a0") and mixture_parameters.containsElementNamed("b0")) {
				return  new Mixture_MultivariateBernoulli (mixture_parameters["a0"], mixture_parameters["b0"]);
			}

			return NULL;
		}
		Prior* gen_prior (Rcpp::List mix_components_prior, Rcpp::List  mix_weight_prior , Rcpp::RObject y) {

			// Q is Fixed lambda or Poisson or Binomial
			// H is gamma

			poisson_gamma_q_param_t * q;
			poisson_gamma_h_param_t * h;

			VERBOSE_ASSERT (mix_components_prior.containsElementNamed("type"), "In gen_prior mix_components_prior does not contain a type field.");
			VERBOSE_ASSERT (mix_weight_prior.containsElementNamed("type"), "In gen_prior mix_weight_prior does not contain a type field.");

			std::string prior_type   = mix_components_prior["type"];
			std::string weight_type  = mix_weight_prior["type"];

			int n = 0;
			if(Rcpp::is<Rcpp::NumericVector>(y) || Rcpp::is<Rcpp::IntegerVector>(y)){
				n = Rcpp::as<arma::vec>(y).size();
			} else  if (Rcpp::is<Rcpp::NumericMatrix>(y) || Rcpp::is<Rcpp::IntegerMatrix>(y)) {
				n = Rcpp::as<arma::mat>(y).n_rows;
			} else {
				VERBOSE_ERROR ("Unsupported type: y variable should be Matrix or Vector.");
			}


			/// ************ Q prior_components_type *******************

			if (mix_components_prior.containsElementNamed("init")
					and mix_components_prior.containsElementNamed("a")
					and mix_components_prior.containsElementNamed("b")) {
				q = new poisson_gamma_q_param_t(
						Rcpp::as<double>(mix_components_prior["init"]),
						Rcpp::as<double>(mix_components_prior["a"]),
						Rcpp::as<double>(mix_components_prior["b"]));
			} else if (mix_components_prior.containsElementNamed("a")
					and mix_components_prior.containsElementNamed("b")) {
					q = new poisson_gamma_q_param_t(
							Rcpp::as<double>(mix_components_prior["a"]),
							Rcpp::as<double>(mix_components_prior["b"]));
			} else if (mix_components_prior.containsElementNamed("Lambda")) {
					q = new poisson_gamma_q_param_t(
							Rcpp::as<double>(mix_components_prior["Lambda"]));
			} else {
				q = new poisson_gamma_q_param_t(n / 10);
			}



			if (mix_weight_prior.containsElementNamed("init")
							and mix_weight_prior.containsElementNamed("a")
							and mix_weight_prior.containsElementNamed("b")) {
				h = new poisson_gamma_h_param_t(
					Rcpp::as<double>(mix_weight_prior["init"]),
					Rcpp::as<double>(mix_weight_prior["a"]),
					Rcpp::as<double>(mix_weight_prior["b"]),
					1);
			} else if (mix_weight_prior.containsElementNamed("a")
					and mix_weight_prior.containsElementNamed("b")) {
			h = new poisson_gamma_h_param_t(
					1,
					Rcpp::as<double>(mix_weight_prior["a"]),
					Rcpp::as<double>(mix_weight_prior["b"]),
					1);
			} else if (mix_weight_prior.containsElementNamed("gamma")) {
				h = new poisson_gamma_h_param_t(
					Rcpp::as<double>(mix_weight_prior["gamma"]));
			} else {
				h = new poisson_gamma_h_param_t(1);
			}



			VERBOSE_ASSERT(q, "Prior error : q is not set.");
			VERBOSE_ASSERT(h, "Prior error : h is not set.");



			Prior* 	prior = new  PriorPoissonGamma (*h,*q);


			return prior;
		}

// INTERNAL FUNCTION - No DOCUMENTATION, Please only use AM_ refixed function.
// [[Rcpp::export]]
Rcpp::List IAM_mcmc_fit (
		Rcpp::RObject                        y                      , /* Not optional */
		Rcpp::List                           mix_kernel_hyperparams , /* Not optional */
		Rcpp::IntegerVector                  initial_clustering     , //  = Rcpp::IntegerVector::create() /* default will be 1for1 */
		Rcpp::List                           mix_components_prior   , //  = Rcpp::List::create()          /* default will be Poisson ()  */
		Rcpp::List                           mix_weight_prior       , //  = Rcpp::List::create()          /* default will be Gamma ()    */
		Rcpp::List                           mcmc_parameters          //  = Rcpp::List::create()          /* (default niter=20000, ….)   */
		) {


	VERBOSE_ASSERT(mcmc_parameters.containsElementNamed("niter"), "mcmc_parameters does not contains niter field." );
	VERBOSE_ASSERT(mcmc_parameters.containsElementNamed("burnin"), "mcmc_parameters does not contains burnin field." );
	VERBOSE_ASSERT(mcmc_parameters.containsElementNamed("thin"), "mcmc_parameters does not contains thin field." );
	VERBOSE_ASSERT(mcmc_parameters.containsElementNamed("verbose"), "mcmc_parameters does not contains verbose field." );
	VERBOSE_ASSERT(mcmc_parameters.containsElementNamed("output"), "mcmc_parameters does not contains output field." );
	VERBOSE_ASSERT(mcmc_parameters.containsElementNamed("file_output"), "mcmc_parameters does not contains file_output field." );

	VERBOSE_LEVEL = mcmc_parameters["verbose"] ;

	VERBOSE_DEBUG ("Debug mode is on.");
	VERBOSE_INFO ("Start mcmc_fit");
	VERBOSE_INFO ("- y = " << y);
	VERBOSE_INFO ("- mix_kernel_hyperparams = " << mix_kernel_hyperparams.size());
	VERBOSE_INFO ("- initial_clustering = " << initial_clustering.size());
	VERBOSE_INFO ("- mix_components_prior = " << mix_components_prior.size());
	VERBOSE_INFO ("- mix_weight_prior = " << mix_weight_prior.size());
	VERBOSE_INFO ("- mcmc_parameters = " << mcmc_parameters.size());
	VERBOSE_INFO ("- Rcpp::is<Rcpp::NumericVector>(y) = " << Rcpp::is<Rcpp::NumericVector>(y));
	VERBOSE_INFO ("- Rcpp::is<Rcpp::IntegerVector>(y) = " << Rcpp::is<Rcpp::IntegerVector>(y));
	VERBOSE_INFO ("- Rcpp::is<Rcpp::NumericMatrix>(y) = " << Rcpp::is<Rcpp::NumericMatrix>(y));
	VERBOSE_INFO ("- Rcpp::is<Rcpp::IntegerMatrix>(y) = " << Rcpp::is<Rcpp::IntegerMatrix>(y));




	Prior*                prior    = gen_prior(mix_components_prior,mix_weight_prior, y);
	Mixture*              mixture  = gen_mix(mix_kernel_hyperparams);


	VERBOSE_ASSERT(mixture, "gen_mix returned NULL");
	VERBOSE_ASSERT(prior, "gen_prior returned NULL");



	if (Rcpp::is<Rcpp::NumericMatrix>(y) || Rcpp::is<Rcpp::IntegerMatrix>(y)) {
			VERBOSE_ASSERT (is_multivariate(mix_kernel_hyperparams), "y argument is a Matrix while the technique is not MultiVariate.") ;

			GibbsResult res =  dynamic_cast<MultivariateMixture*>(mixture)->fit(Rcpp::as<arma::mat>(y) , initial_clustering, prior ,
					mcmc_parameters["niter"] ,mcmc_parameters["burnin"] ,mcmc_parameters["thin"] ,mcmc_parameters["verbose"] );

			return Rcpp::List::create(
						Rcpp::Named("U_post")      = res.U,
						Rcpp::Named("ci_post")     = res.ci,
						Rcpp::Named("S_post")      = res.S,
						Rcpp::Named("M_post")      = res.M,
						Rcpp::Named("K_post")      = res.K);
		}  else if(Rcpp::is<Rcpp::NumericVector>(y) || Rcpp::is<Rcpp::IntegerVector>(y)){
			VERBOSE_ASSERT (is_univariate(mix_kernel_hyperparams), "y argument is a Vector while the technique is not Univariate.") ;
			GibbsResult res =  dynamic_cast<UnivariateMixture*>(mixture)->fit(Rcpp::as<arma::vec>(y) , initial_clustering, prior ,
					mcmc_parameters["niter"] ,mcmc_parameters["burnin"] ,mcmc_parameters["thin"] ,mcmc_parameters["verbose"] );

			return Rcpp::List::create(
						Rcpp::Named("U_post")      = res.U,
						Rcpp::Named("ci_post")     = res.ci,
						Rcpp::Named("S_post")      = res.S,
						Rcpp::Named("M_post")      = res.M,
						Rcpp::Named("K_post")      = res.K);

		} else {
		VERBOSE_ERROR("The parameter y must be a Matrix or a Vector.");
	}
	return Rcpp::List::create(Rcpp::Named("Error") = "Unexpected error."  );

}
