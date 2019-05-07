/*
 * PriorDirac.hpp
 *
 *  Created on: Apr 12, 2019
 *      Author: toky
 */

#ifndef ANTMAN_SRC_PRIORDIRAC_HPP_
#define ANTMAN_SRC_PRIORDIRAC_HPP_


#include <cassert>
#include <RcppArmadillo.h> // [[Rcpp::depends(RcppArmadillo)]]
#include "Prior.hpp"

// --------------------------------------------------------------------------------------------------------------------

struct dirac_q_param_t {
	double Mstar ;
	dirac_q_param_t (double Mstar) : Mstar (Mstar) {}
};

struct dirac_h_param_t {
	double gamma ;
	dirac_h_param_t (double gamma) : gamma (gamma) {}
};

// --------------------------------------------------------------------------------------------------------------------

class PriorDirac : public TypedPrior < dirac_h_param_t, dirac_q_param_t> {

private:

	double ah,bh; //hyper-prior parameters for h
	double lsd, lsd_g; //this is the standard deviation of the MH algorithm to update gamma.



public:
	PriorDirac (double inith, double initq, double ah, double bh , double lsd) :
		TypedPrior<dirac_h_param_t, dirac_q_param_t>(inith,initq,false),
		ah(ah), bh(bh) ,
	 lsd(lsd) , lsd_g(1)   {
			// TODO[CHECK ME] : requirements ? what should I check.
		    // TODO [Unfinished business] If I don't have init, generate
	}

	double get_gamma() {return this->h_param.gamma;};

	void init_q_param () {
	}
	void init_h_param () {
		this->h_param.gamma = R::rgamma(ah,bh);
		VERBOSE_ASSERT(this->h_param.gamma > 0, "Please provide gamma, R::rgamma(ah,bh) returned 0.");
	}



	void update_q_param (const  double U, const  int K) {}

	inline double log_full_gamma( const double Loc_gamma, const int K , const  std::vector<int> & nj,const   int Mstar,const  double U_current ,const  double ag,const  double bg){

		double out= - ( Loc_gamma  * Mstar ) * std::log (1 + U_current) ;



		for(int j=0;j<K;j++){
			out+=std::lgamma(Loc_gamma+ (double) nj[j])-std::lgamma(Loc_gamma);
		}
		/// When the prior is a gamma

		out+=(ag-1)*std::log(Loc_gamma)-bg*Loc_gamma;

		return(out);
	}




	void update_h_param (const  double U, const  int K, const std::vector<int> &nj ) {

			const double vecchio = this->h_param.gamma;
			const double lmedia = std::log(vecchio);

			//Propose a new value
			const double lnuovo=R::rnorm(lmedia,lsd);
			const double nuovo=std::exp(lnuovo);

			double ln_acp = log_full_gamma(nuovo , K , nj,  this->q_param.Mstar, U , ah, bh ) - lmedia;

			ln_acp= ln_acp - (log_full_gamma(vecchio , K , nj,   this->q_param.Mstar, U , ah, bh) - lnuovo);

			const double lnu=std::log(R::runif(0.0,1.0));


			this->h_param.gamma = lnu<ln_acp ? nuovo : vecchio;
			VERBOSE_DEBUG( " (inside Dirac) LSD is " << lsd);
			lsd = update_lsd (  lsd,  ln_acp,  lsd_g++) ;

	}


	int init_M_na(const int K) {
		int M_na = this->q_param.Mstar - K ;
		VERBOSE_ASSERT(M_na >= 0, "Please provide initial clustering with K <= Mstar");
		return M_na;
	}
	int update_M_na(const double U ,  const int K) {
		int M_na = this->q_param.Mstar - K ;
		VERBOSE_ASSERT(M_na >= 0, "Internal Error, K > Mstar.");
		return M_na;

	}

};
// --------------------------------------------------------------------------------------------------------------------




#endif /* ANTMAN_SRC_PRIORDIRAC_HPP_ */
