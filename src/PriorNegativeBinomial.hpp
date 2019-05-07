/*
 * MixturePrior.hpp
 *
 *  Created on: Apr 8, 2019
 *      Author: toky
 */

#ifndef ANTMAN_SRC_PRIORNEGATIVEBINOMIAL_HPP_
#define ANTMAN_SRC_PRIORNEGATIVEBINOMIAL_HPP_


#include <cassert>
#include <RcppArmadillo.h> // [[Rcpp::depends(RcppArmadillo)]]
#include "Prior.hpp"
#include "utils.hpp"

// --------------------------------------------------------------------------------------------------------------------

struct negative_binomial_q_param_t {
	double R_M, P_M;
	negative_binomial_q_param_t ()  : R_M (0), P_M (0)  {}
	negative_binomial_q_param_t (double R_M, double P_M) : R_M (R_M), P_M (P_M) {}
};

struct negative_binomial_h_param_t {
	double gamma ;
	negative_binomial_h_param_t ()  : gamma(0)  {}
	negative_binomial_h_param_t (double gamma) : gamma (gamma) {}
};

// --------------------------------------------------------------------------------------------------------------------

class PriorNegativeBinomial : public TypedPrior < negative_binomial_h_param_t, negative_binomial_q_param_t> {

private:

	double a_R, b_R; //hyper-prior parameters for q
	double a_P, b_P; //hyper-prior parameters for q

	double LSDR_M, LSDP_M;
	double LSDR_M_g, LSDP_M_g;

	double a,b; //hyper-prior parameters for h
	double lsd,lsd_g; //this is the standard deviation of the MH algorithm to update gamma.

	PriorNegativeBinomial () :
		TypedPrior<negative_binomial_h_param_t, negative_binomial_q_param_t>(),
		a_R(0.0), b_R(0.0) , a_P(0.0), b_P(0.0),
		LSDR_M(0.0) , LSDP_M(0.0) ,
		LSDR_M_g(0.0) , LSDP_M_g(1.0) ,
		 a(0.0), b(0.0), lsd(0.0), lsd_g(1.0) {};

public:
	PriorNegativeBinomial (double initgamma,
			double R_M, double P_M ,
			double LSDR_M, double LSDP_M, double a_R, double b_R, double a_P, double b_P , double a, double b, double lsd) :
		TypedPrior<negative_binomial_h_param_t, negative_binomial_q_param_t>(negative_binomial_h_param_t(initgamma), negative_binomial_q_param_t(R_M,P_M),false)
		, a_R(a_R), b_R(b_R) , a_P(a_P), b_P(b_P) ,
		LSDR_M(LSDR_M) , LSDP_M(LSDP_M),
		LSDR_M_g(0.0) , LSDP_M_g(0.0) ,
		 a(a), b(a), lsd(lsd), lsd_g(0.0)  {
			assert ((P_M <=1) && (P_M >=0));
	}

	void init_q_param () {

	}
	void init_h_param () {
		this->h_param.gamma = R::rgamma(a,b);
		VERBOSE_ASSERT(this->h_param.gamma > 0, "Please provide gamma, R::rgamma(ah,bh) returned 0.");
	}


	inline double compute_lphi (double U_current, double h_param) {
		double lphi_u =  - h_param * std::log( 1 + U_current ) ;
		return lphi_u;

	}

	inline double compute_lkappa ( double U_current , double h_param, double K , int nij ,  double P_M , double R_M) {
		double log_kappa = - ( nij + h_param ) * std::log (U_current + 1) +  std::lgamma( h_param + (double) nij) - std::lgamma(h_param) ;

		return log_kappa;
	}

	inline double compute_lPsi ( double U_current , double h_param, double K ,  double P_M , double R_M) {

		// h is a gamma density here
		double phi_u = std::exp( compute_lphi(U_current, h_param) );

		double lfirst = std::lgamma(R_M + K - 1) - std::lgamma (R_M) + (K - 1) * std::log (P_M) + R_M * std::log (1 - P_M) ;


		double log_num   =  std::log ( phi_u * (R_M -1 ) + K ) ;
		double log_denum = ( K + R_M) * std::log (1 - phi_u * P_M ) ;

		double lpsi = lfirst + log_num - log_denum;

		return lpsi;
	}

	inline double log_full_EPPF(  double h_param, double K , const std::vector<int> & nj, double U_current ,
			double P_M , double R_M ){
		// TODO[CHECK ME] : Take care this is proportionnal to U
		double lPsi = compute_lPsi ( U_current , h_param, K ,  P_M , R_M) ;

		double log_kappa_sum = 0 ;
		for(int j=0;j<K;j++){
			double log_kappa = compute_lkappa (U_current , h_param, K , nj[j] ,  P_M , R_M);
			log_kappa_sum += log_kappa;
		}

		double out = lPsi + log_kappa_sum;

		return out;
	}


void update_q_param (const  double U, const  int K) {

		const double h_param = this->h_param.gamma;
		const double R_M     = this->q_param.R_M;
		const double P_M     = this->q_param.P_M;


		// Metropolis-Hasting for R_M

		double R_vecchio = R_M;
		double R_lmedia = std::log(R_vecchio);

		//Propose a new value
		double R_lnuovo=R::rnorm(R_lmedia,LSDR_M);
		double R_nuovo=std::exp(R_lnuovo);

		double log_full_r_m_new =  compute_lPsi ( U ,  h_param,  K ,   P_M ,  R_nuovo) + (a_R-1)*std::log(R_nuovo)-b_R*R_nuovo ;
		double log_full_r_m_vec =  compute_lPsi ( U ,  h_param,  K ,   P_M ,  R_vecchio) + (a_R-1)*std::log(R_vecchio)-b_R*R_vecchio;

		double R_ln_acp = (log_full_r_m_new - R_lmedia) - (log_full_r_m_vec - R_lnuovo);

		double R_lnu=std::log(R::runif(0.0,1.0));

		this->q_param.R_M = R_lnu < R_ln_acp ? R_nuovo : R_vecchio;

		LSDR_M = update_lsd (  LSDR_M,  R_ln_acp,  LSDR_M_g++) ;


		// Metropolis-Hasting for P_M
		// TODO[CHECK ME] : (Raffa should check this ... at some point .. in time ... we are done ...)

		double P_vecchio = P_M;
		double P_lmedia = std::log(P_vecchio) - std::log(1 - P_vecchio);

		//Propose a new value
		double P_lnuovo=R::rnorm(P_lmedia,LSDP_M);
		double P_nuovo=std::exp(P_lnuovo) / (1 + std::exp (P_lnuovo));

		double log_full_p_m_new =  compute_lPsi ( U ,  h_param,  K ,   P_M ,  R_nuovo) + (a_R-1)*std::log(P_nuovo)+(b_R - 1)* std::log(1 - P_nuovo) ;
		double log_full_p_m_vec =  compute_lPsi ( U ,  h_param,  K ,   P_M ,  R_vecchio) + (a_R-1)*std::log(P_vecchio)-(b_R - 1 )*std::log(1 -P_vecchio);

		double P_ln_acp = (log_full_p_m_new - P_lmedia - std::log(1 - P_vecchio ) ) - (log_full_p_m_vec - P_lnuovo - std::log(1 - P_nuovo ));

		double P_lnu=std::log(R::runif(0.0,1.0));

		this->q_param.P_M = P_lnu < P_ln_acp ? P_nuovo : P_vecchio;
		LSDP_M = update_lsd (  LSDP_M,  P_ln_acp,  LSDP_M_g++) ;

		return;
	}

void update_h_param (const double U ,  const int K , const  std::vector<int> & nj ) {

		const double R_M     = this->q_param.R_M;
		const double P_M     = this->q_param.P_M;
		const double vecchio = this->h_param.gamma;
		const double lmedia = std::log(vecchio);

		//Propose a new value
		const double lnuovo=R::rnorm(lmedia,lsd);
		const double nuovo=std::exp(lnuovo);

		const double log_full_gamma_new = log_full_EPPF (nuovo , K , nj,  U , P_M, R_M ) + (a-1)*std::log(nuovo)-b*nuovo;
		const double log_full_gamma_vec = log_full_EPPF (vecchio , K , nj,  U , P_M, R_M ) + (a-1)*std::log(vecchio)-b*vecchio;

		const double ln_acp = (log_full_gamma_new - lmedia) - (log_full_gamma_vec - lnuovo);

		const double lnu=std::log(R::runif(0.0,1.0));

		this->h_param = lnu < ln_acp ? nuovo : vecchio;

		lsd = update_lsd (  lsd,  ln_acp,  lsd_g++) ;



	}
int init_M_na() {
		const double R_M     = this->q_param.R_M;
		const double P_M     = this->q_param.P_M;
		// TODO[CHECK ME]: According to Rafaelle what we call P is 1-1 in R. and what we call R_M is Size in R.
		return R::rnbinom(R_M, 1-P_M);
	}

int update_M_na(const double U ,  const int K) {

		const double R_M     = this->q_param.R_M;
		const double P_M     = this->q_param.P_M;

		int M_na;

		const double phi_u  = 1 / std::pow(1 + U, this->h_param.gamma);
		const double lphi_u = - this->h_param.gamma*log(1+U);
		const double up     = (R_M + K) * phi_u * (P_M) ;
		const double down   =  R_M * phi_u *P_M + K ;
		const double peso   =  up / down;

		const double unif=R::runif(0.0,1.0);

		if(unif < peso){
			M_na=R::rnbinom(R_M + K, phi_u * P_M) + 1;
		} else {
			M_na=R::rnbinom(R_M - 1 + K, phi_u * P_M) ;
		}
		return M_na;

	}


	double get_gamma() {return this->h_param.gamma;};

};
// --------------------------------------------------------------------------------------------------------------------

#endif /* ANTMAN_SRC_PRIORNEGATIVEBINOMIAL_HPP_ */
