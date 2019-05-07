/*
 * Mixture.hpp
 *
 *  Created on: 14 Apr 2019
 *      Author: toky
 */

#ifndef ANTMAN_SRC_MIXTURE_HPP_
#define ANTMAN_SRC_MIXTURE_HPP_

#include "utils.hpp"
#include "Prior.hpp"
#include "GibbsResult.hpp"


class allocation_result {

	cluster_indices_t          _ci;
	std::vector<int>           _nj;
	Rcpp::NumericVector         _S;
public :
	allocation_result (cluster_indices_t & ci , std::vector<int>  & nj, Rcpp::NumericVector & S) : _ci (ci) , _nj(nj), _S(S) {} ;

	inline const Rcpp::NumericVector  &        S          () const {return _S;}
	inline const std::vector<int>     &        nj         () const {return _nj;}
	inline const cluster_indices_t    &        ci         () const {return _ci;}
};

class Mixture {
public :
	virtual    ~Mixture() {};
};


template<typename InputType>
class TypedMixture : public Mixture {
public :
	virtual    ~TypedMixture() {};
protected:
	typedef InputType input_t ;

	virtual void init_tau (const InputType & y, const int M) = 0;

	virtual cluster_indices_t  up_ci(const  InputType & y,
			const long M,
			const Rcpp::NumericVector & S) = 0;

	 virtual allocation_result up_allocated_nonallocated (
			const int K ,
			const int M ,
			const cluster_indices_t & ci ,
			const cluster_indices_t & ci_star  ,
			const double gamma,
			const double U,
			const  InputType & y ) = 0;


public:

	GibbsResult fit(InputType y,
			cluster_indices_t initial_clustering,
			Prior * prior,
			const int niter,
			const int burnin ,
			const int thin,
			const int verbose ) {


		const int n = y.n_rows;

		/*******************************************/
		/*******       State Variables            **/
		/*******************************************/

		double               U_current      = 0;
		cluster_indices_t    ci_current     = initial_clustering;

		// Then the number of allocated jumps (i.e. number of cluster) is
		cluster_indices_t ci_star = arma::unique(ci_current);
		int K=ci_star.size();


		VERBOSE_DEBUG("prior->init();");
		prior->init();
		VERBOSE_DEBUG("Done");

		int M_na= prior->init_M_na();
		int M=K+M_na;
		VERBOSE_DEBUG("this->init_tau (y, M);");
		this->init_tau (y, M);
		VERBOSE_DEBUG("Done");

		Rcpp::NumericVector  S_current=Rcpp::NumericVector(M);

		// TODO : S current initialized with gamma_current !!! should not right ??
		for(int it=0;it<M;it++){
			S_current[it] =R::rgamma(2.0,1.0); // replace gamma current by 2.0 for now
		}



		/**************************************/
		/*******       initialize the output **/
		/**************************************/


		GibbsResult result;

		result.U.resize(niter);
		result.ci.resize(niter);
		result.S.resize(niter);
		result.M.resize(niter);
		result.K.resize(niter);

		/**************************************/
		/******* Start Gibbs                 **/
		/**************************************/

		Rcpp::Rcout <<"Let's start the Gibbs!"<<"\n";

		double total_gibbs           = 0;
		auto   start_gibbs           = std::chrono::system_clock::now();

		int iter=0;

		while (iter < (niter+burnin)) {

			double total_iter           = 0;
			double total_u              = 0;
			double total_ci             = 0;
			double total_mna            = 0;
			double total_alloc          = 0;
			auto start_iter = std::chrono::system_clock::now();

			for(int thi=0;thi<thin;thi++){

				//Following the scheme in Figure 1 of the paper
				//I first update the latent U

				auto start_u = std::chrono::system_clock::now();
				VERBOSE_DEBUG("S_current = " << S_current <<"\n");
				U_current=R::rgamma(n,1.0/Rcpp::sum(S_current));
				VERBOSE_DEBUG("U_current = " << U_current <<"\n");
				auto end_u = std::chrono::system_clock::now();
				auto elapsed_u = end_u - start_u;
				total_u += elapsed_u.count() / 1000000.0;


				// Update CI and CI*
				VERBOSE_DEBUG("Call up_ci\n");
				auto start_ci = std::chrono::system_clock::now();				
				if(iter>0){//I need this if i want that my inizialization for ci works
				ci_current = this->up_ci(y, M, S_current); // parametricPrior = k_x,,X   Tau = Beta_current, z_current
				}
				// TODO: This is computed twice, could be avoided.
				cluster_indices_t ci_star = arma::unique(ci_current);
				K=ci_star.size();
				auto end_ci = std::chrono::system_clock::now();
				auto elapsed_ci = end_ci - start_ci;
				total_ci += elapsed_ci.count() / 1000000.0;
				VERBOSE_DEBUG("End up_ci\n");


				auto start_mna = std::chrono::system_clock::now();
				M_na = prior->update_M_na(U_current, K);
				M=K+M_na;
				auto end_mna = std::chrono::system_clock::now();
				auto elapsed_mna = end_mna - start_mna;
				total_mna += elapsed_mna.count() / 1000000.0;


				VERBOSE_DEBUG ( "K= " << K << "M= " << M  << std::endl);
				VERBOSE_DEBUG ( "ci_star="<< ci_star << std::endl);
				VERBOSE_DEBUG ( "ci_current="<< ci_current << std::endl);
				VERBOSE_DEBUG ( "gamma_current="<< prior->get_gamma() << std::endl);
				VERBOSE_DEBUG ( "U_current="<< U_current << std::endl);


				// Compute Allocation (ci_reorder, nj, Beta and S (allocSide))

				VERBOSE_DEBUG("Call up_allocated_nonallocated\n");
				auto start_alloc = std::chrono::system_clock::now();
				auto up_allocated_res = this->up_allocated_nonallocated (  K , M , ci_current ,  ci_star  , prior->get_gamma(),  U_current, y );
				const std::vector <int> nj = up_allocated_res.nj();
				assert(K == nj.size());
				ci_current   = up_allocated_res.ci();
				S_current = up_allocated_res.S();
				auto end_alloc = std::chrono::system_clock::now();
				auto elapsed_alloc = end_alloc - start_alloc;
				total_alloc += elapsed_alloc.count() / 1000000.0;
				VERBOSE_DEBUG("End up_allocated_nonallocated\n");

				prior->update(U_current, K, nj);

			}

			auto end_iter = std::chrono::system_clock::now();
			auto elapsed_iter = end_iter - start_iter;
			total_iter = elapsed_iter.count()  / 1000000.0;

			if(verbose!=0){
				if((iter% verbose)==0){

					auto end_gibbs             = std::chrono::system_clock::now();
					auto elapsed_gibbs         = end_gibbs - start_gibbs;
					     start_gibbs           = std::chrono::system_clock::now();
					total_gibbs                = elapsed_gibbs.count() / 1000000.0;

					Rcpp::Rcout <<"iter="<<iter<<" K="<<K<<" M_na="<<M_na<<" M="<<M<<
							" u=" <<total_u <<
							"ms ci=" <<total_ci <<
							"ms mna=" <<total_mna <<
							"ms alloc=" <<total_alloc <<
							"ms total_iter=" <<total_iter <<
							"ms total_gibbs=" <<total_gibbs<< "ms" << std::endl;

				}
			}


			// Save output after the burn-in and taking into account
			// thinning
			if( (iter>=burnin)){

				result.U[iter-burnin]=U_current;

				// NO NON NO DEVO STARE MOLTO ATTENTO ci_out[iter-burnin]=ci_current;



				result.ci[iter-burnin]=Rcpp::IntegerVector(n);
				std::copy(ci_current.begin(),ci_current.end(),result.ci[iter-burnin].begin());
				result.S[iter-burnin]=S_current;
				result.M[iter-burnin]=M;
				result.K[iter-burnin]=K;
			}



			iter +=1;

		}//I close the while

		Rcpp::Rcout <<"End of Iterations." << std::endl;
		return result;

	}

};

class UnivariateMixture : public TypedMixture<arma::vec> {};
class MultivariateMixture : public TypedMixture<arma::mat> {};

#endif /* ANTMAN_SRC_MIXTURE_HPP_ */
