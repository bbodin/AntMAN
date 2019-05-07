/*
 * PoissonPrior.hpp
 *
 *  Created on: Apr 12, 2019
 *      Author: toky
 */

#ifndef ANTMAN_SRC_PRIORPOISSON_HPP_
#define ANTMAN_SRC_PRIORPOISSON_HPP_


#include <cassert>
#include <RcppArmadillo.h> // [[Rcpp::depends(RcppArmadillo)]]
#include "Prior.hpp"

// --------------------------------------------------------------------------------------------------------------------


template<typename q_param_t>
class gamma_h_param_t {
public:
	bool gamma_is_fixed;
	double gamma;
	const double a,b; // hyper-prior parameters for h
	double lsd, lsd_g;       // this is the standard deviation of the MH algorithm to update gamma.
	gamma_h_param_t (double gamma, double a, double b, double lsd) : gamma_is_fixed (false), gamma (gamma) , a(a), b(b), lsd(lsd), lsd_g(1) {}
	gamma_h_param_t (              double a, double b, double lsd) : gamma_is_fixed (false), gamma (R::rgamma(a,b)) , a(a), b(b), lsd(lsd), lsd_g(1) {}
	gamma_h_param_t (double gamma) : gamma_is_fixed (true), gamma (gamma), a(0), b(0), lsd(0), lsd_g(1)  {}

	inline double log_full_gamma( const double Loc_gamma, const int K ,const  std::vector<int> & nj,const   double Lambda_current,const  double U_current ,const  double ag,const  double bg){

		double out=0;
		double up1g=std::pow(1+U_current,Loc_gamma);

		out+=std::log(Lambda_current/up1g+K)+Lambda_current/up1g-K*std::log(up1g);
		for(int j=0;j<K;j++){
			out+=std::lgamma(Loc_gamma+ (double) nj[j])-std::lgamma(Loc_gamma);
		}
		/// When the prior is a gamma

		out+=(ag-1)*std::log(Loc_gamma)-bg*Loc_gamma;

		return(out);
	}

	void update (const  double U, const  int K, const std::vector<int> &nj , const q_param_t& q_param) {
		if (this->gamma_is_fixed) return;

		const double vecchio = this->gamma;
		const double lmedia = std::log(vecchio);

			//Propose a new value
			const double lnuovo=R::rnorm(lmedia,lsd);
			const double nuovo=std::exp(lnuovo);


	//		const double log_full_gamma_new = log_full_EPPF (nuovo , K , nj,  U , this->q_param.lambda ) + (ah-1)*std::log(nuovo)-bh*nuovo;
	//		const double log_full_gamma_vec = log_full_EPPF (vecchio , K , nj,  U , this->q_param.lambda ) + (ah-1)*std::log(vecchio)-bh*vecchio;
	//		const double ln_acp = (log_full_gamma_new - lmedia) - (log_full_gamma_vec - lnuovo);


			double ln_acp = log_ful1l_gamma(nuovo , K , nj,  q_param.lambda, U , a, b ) - lmedia;

			ln_acp= ln_acp - (log_full_gamma(vecchio , K , nj,   q_param.lambda, U , a, b) - lnuovo);

			const double lnu=std::log(R::runif(0.0,1.0));

			this->gamma = lnu<ln_acp ? nuovo : vecchio;

			lsd = update_lsd (  lsd,  ln_acp,  lsd_g++) ;

	}

};

class poisson_gamma_q_param_t {
public:
	bool lambda_is_fixed;
	double lambda ;
	const double a,b; // hyper-prior parameters for q
	poisson_gamma_q_param_t (double lambda, double a, double b) : lambda_is_fixed(false), lambda (lambda) ,a(a),b(b) {}
	poisson_gamma_q_param_t (               double a, double b) : lambda_is_fixed(false), lambda (1.0)    ,a(a),b(b) {}
	poisson_gamma_q_param_t (double lambda)                     :  lambda_is_fixed(true), lambda (lambda) ,a(0),b(0) {}

	void update (const  double U, const  int K, const gamma_h_param_t <poisson_gamma_q_param_t>& h_param) {
		if (lambda_is_fixed) return;
		double lphi_u=-h_param.gamma*std::log(1+U);
		double lpeso=lphi_u-std::log(1+b);
		double lunif=std::log(R::runif(0.0,1.1));
		double astar=K+a ;// See the notation of Point 4 in 10.1
		double rate=(1-std::exp(lphi_u)+b);
		this->lambda = lunif<lpeso ? R::rgamma(astar+1,1.0/rate) : R::rgamma(astar,1.0/rate);
	}


};

typedef gamma_h_param_t<poisson_gamma_q_param_t> poisson_gamma_h_param_t ;



class PriorPoissonGamma : public TypedPrior <poisson_gamma_h_param_t, poisson_gamma_q_param_t> {
public:
	PriorPoissonGamma (poisson_gamma_h_param_t gamma_h_param, poisson_gamma_q_param_t poisson_q_param) :
		TypedPrior <poisson_gamma_h_param_t, poisson_gamma_q_param_t> (gamma_h_param,poisson_q_param)  {
		}

	int init_M_na () {
		return R::rpois(this->q_param.lambda);
	}

	int update_M_na (const double U ,  const int K) {

		int M_na;

		const double phi_u=-this->h_param.gamma*log(1+U);
		const double Lambda_u=exp(std::log(this->q_param.lambda)+phi_u);

		const double unif=R::runif(0.0,1.0);

		if(unif<(Lambda_u/(Lambda_u+K))){
			M_na=R::rpois(Lambda_u)+1;

		} else {
			M_na=R::rpois(Lambda_u);
		}
		return M_na;

	}

};

#endif /* ANTMAN_SRC_PRIORPOISSON_HPP_ */
