/*
 * MixturePrior.cpp
 *
 *  Created on: Apr 8, 2019
 *      Author: toky
 */

#include "PriorNegativeBinomial.hpp"

#include <cmath>
#include <vector>


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
		// TODO : Take care this is proportionnal to U
		double lPsi = compute_lPsi ( U_current , h_param, K ,  P_M , R_M) ;

		double log_kappa_sum = 0 ;
		for(int j=0;j<K;j++){
			double log_kappa = compute_lkappa (U_current , h_param, K , nj[j] ,  P_M , R_M);
			log_kappa_sum += log_kappa;
		}

		double out = lPsi + log_kappa_sum;

		return out;
	}


void PriorNegativeBinomial::update_q_param (const  double U, const  int K) {

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



		// Metropolis-Hasting for P_M
		// TODO : (Raffa should check this ... at some point .. in time ... we are done ...)

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

		return;
	}

// TODO : does nj can be called an histogram ?
void PriorNegativeBinomial::update_h_param (const double U ,  const int K , const  std::vector<int> & nj ) {

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

	}
int PriorNegativeBinomial::init_M_na() {
		const double R_M     = this->q_param.R_M;
		const double P_M     = this->q_param.P_M;
		// TODO: According to Rafaelle what we call P is 1-1 in R. and what we call R_M is Size in R.
		return R::rnbinom(R_M, 1-P_M);
	}

int PriorNegativeBinomial::update_M_na(const double U ,  const int K) {

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

