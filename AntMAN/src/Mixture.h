/*
 * Mixture.hpp
 *
 *  Created on: 14 Apr 2019
 */

#ifndef ANTMAN_SRC_MIXTURE_H_
#define ANTMAN_SRC_MIXTURE_H_
#include <map>
#include <iomanip>

#include "GibbsResult.h"
#include "Prior.h"
#include "utils.h"


class allocation_result {

	cluster_indices_t          _ci;
	std::vector<int>           _nj;
	arma::vec                   _S;
public :
	allocation_result (const cluster_indices_t & ci , const  std::vector<int>  & nj, const arma::vec & S) : _ci (ci) , _nj(nj), _S(S) {} ;

	inline const arma::vec            &        S          () const {return _S;}
	inline const std::vector<int>     &        nj         () const {return _nj;}
	inline const cluster_indices_t    &        ci         () const {return _ci;}
};

class Mixture {
public :
	virtual             ~Mixture() {};
#ifdef HAS_RCPP
	virtual Rcpp::List  get_tau () = 0;
#else
	virtual std::string  get_tau() = 0;
#endif
private :
	bool _parallel;
protected :
	void set_parallel(bool parallel) {this->_parallel = parallel;};
	bool get_parallel() {return this->_parallel;};
};


template<typename InputType>
class TypedMixture : public Mixture {
public :
	virtual    ~TypedMixture() {};
protected:
	typedef InputType input_t ;

	virtual void        init_tau (const InputType & y, const int M) = 0;


	virtual cluster_indices_t  up_ci(const  InputType & y,
			const long M,
			const arma::vec & S) = 0;

	virtual allocation_result up_allocated_nonallocated (
			const int K ,
			const int M ,
			const cluster_indices_t & ci ,
			const cluster_indices_t & ci_star  ,
			const double gamma,
			const double U,
			const  InputType & y ) = 0;





public:

	void fit(InputType y,
			cluster_indices_t initial_clustering,
			bool fixed_clustering,
			Prior * prior,
			const unsigned long  niter,
			const unsigned long  burnin ,
			const unsigned long  thin,
			bool parallel,
			GibbsResult * results) {

		VERBOSE_ASSERT(niter > burnin, "Please use a total iteration number greater then burnin.");
		VERBOSE_ASSERT(thin  > 0     , "Please make sure to have thin > 0.");


		this->set_parallel(parallel);

		const int n = y.n_rows;

		/*******************************************/
		/*******       State Variables            **/
		/*******************************************/

		double               U_current      = 0;
		cluster_indices_t    ci_current     = initial_clustering;

		// Then the number of allocated jumps (i.e. number of cluster) is
		cluster_indices_t ci_star = arma::unique(ci_current);
		unsigned int K=ci_star.size();


		int M_na= prior->init_M_na(K);
		int M=K+M_na;
		VERBOSE_DEBUG("this->init_tau (y, M);");
		this->init_tau (y, M);
		VERBOSE_DEBUG("Done");

		arma::vec  S_current (M);

		// TODO[CHECK ME] : S current initialized with gamma_current !!! should not right ??
		for(int it=0;it<M;it++){
			S_current[it] =am_rgamma(2.0,1.0); // replace gamma current by 2.0 for now
		}





		/**************************************/
		/******* Start Gibbs                 **/
		/**************************************/

		VERBOSE_INFO ( "Let's start the Gibbs!");

		unsigned long total_saved          = 0;
		double total_iter           = 0;
		double total_u              = 0;
		double total_ci             = 0;
		double total_mna            = 0;
		double total_alloc          = 0;
		double total_gibbs          = 0;
		auto   start_gibbs          = std::chrono::system_clock::now();


		for (unsigned int iter = 0 ; iter < niter ; iter++)  {

			total_iter           = 0;
			total_u              = 0;
			total_ci             = 0;
			total_mna            = 0;
			total_alloc          = 0;

			auto start_iter = std::chrono::system_clock::now();


			//Following the scheme in Figure 1 of the paper
			//I first update the latent U

			auto start_u = std::chrono::system_clock::now();
			VERBOSE_EXTRA("S_current = " << S_current <<"\n");
			U_current=am_rgamma(n,1.0/arma::sum(S_current));
			VERBOSE_EXTRA("U_current = " << U_current <<"\n");
			auto end_u = std::chrono::system_clock::now();
			auto elapsed_u = end_u - start_u;
			total_u += elapsed_u.count() / 1000000.0;


			// Update CI and CI*
			VERBOSE_DEBUG("Call up_ci\n");
			auto start_ci = std::chrono::system_clock::now();
			if ((iter > 0) and (not fixed_clustering)) {//I need this if i want that my inizialization for ci works
				ci_current = this->up_ci(y, M, S_current); // parametricPrior = k_x,,X   Tau = Beta_current, z_current
			}
			// TODO[OPTIMIZE ME]: This is computed twice, could be avoided.
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
			VERBOSE_EXTRA ( "ci_star="<< ci_star << std::endl);
			VERBOSE_EXTRA ( "ci_current="<< ci_current << std::endl);
			VERBOSE_DEBUG ( "gamma_current="<< prior->get_gamma() << std::endl);
			VERBOSE_EXTRA ( "U_current="<< U_current << std::endl);


			// Compute Allocation (ci_reorder, nj, Beta and S (allocSide))

			VERBOSE_DEBUG("Call up_allocated_nonallocated\n");
			auto start_alloc = std::chrono::system_clock::now();
			auto up_allocated_res = this->up_allocated_nonallocated (  K , M , ci_current ,  ci_star  , prior->get_gamma(),  U_current, y );
			const std::vector <int> nj = up_allocated_res.nj();

			ci_current   = up_allocated_res.ci();

			S_current = up_allocated_res.S();
			auto end_alloc = std::chrono::system_clock::now();
			auto elapsed_alloc = end_alloc - start_alloc;
			total_alloc += elapsed_alloc.count() / 1000000.0;
			VERBOSE_DEBUG("End up_allocated_nonallocated\n");

			prior->update(U_current, K, nj);

			VERBOSE_DEBUG("prior->update(U_current, K, nj) is done\n");

			auto end_iter = std::chrono::system_clock::now();
			auto elapsed_iter = end_iter - start_iter;
			total_iter = elapsed_iter.count()  / 1000000.0;
			VERBOSE_DEBUG("total_iter = " << total_iter << "ms");

			const unsigned long verbose_slice = niter / std::min((long unsigned int)100,niter);
			const unsigned long total_to_save = ((niter - burnin) / thin) + ((((niter - burnin) % thin) == 0)?0:1) ;

			// Save output after the burn-in and taking into account
			// thinning
			if( (iter >= burnin) and ((iter - burnin) % thin == 0) ) {
				total_saved++;
				results->log_output (ci_current,  S_current,  M,  K,  M_na, this , prior) ;
				VERBOSE_ASSERT(total_to_save >= total_saved, "Raffaele was right.");
				VERBOSE_DEBUG("results->log_output() is done");
			} else {
				VERBOSE_DEBUG("results->log_output() is skiped");
			}

			// Logging

			VERBOSE_DEBUG("verbose_slice = " << verbose_slice);
			if((((iter % (verbose_slice))==0) or  ((iter + 1) == niter))) {
				VERBOSE_DEBUG("Start the logging");
				auto end_gibbs             = std::chrono::system_clock::now();
				auto elapsed_gibbs         = end_gibbs - start_gibbs;
				start_gibbs           = std::chrono::system_clock::now();
				total_gibbs               += elapsed_gibbs.count() / 1000000.0;

				VERBOSE_LOG("[" << std::setw(3) << 100 * iter / (niter - 1) << "%]" <<
						" iter=["<< iter + 1 << "/" << niter << "]" <<
						" saved=["<< total_saved << "/" << total_to_save << "]" <<
						" K="<<K<<
						//" M_na="<<M_na<<
						" M="<<M<<
						" ci=" <<total_ci << "ms" <<
						//"ms mna=" <<total_mna << "ms" <<
						" alloc=" <<total_alloc << "ms" <<
						" iter=" <<total_iter << "ms" <<
						" total_gibbs=" <<total_gibbs<< "ms" );
			} else {
				VERBOSE_DEBUG("Skip the logging");
			}

		}//I close the while

		VERBOSE_LOG("End of Iterations." );


	}

};

class UnivariateMixture : public TypedMixture<arma::vec> {};
class MultivariateMixture : public TypedMixture<arma::mat> {};

#endif /* ANTMAN_SRC_MIXTURE_H_ */
