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

	double a,b; //hyper-prior parameters for h
	double lsd; //this is the standard deviation of the MH algorithm to update gamma.

	PriorNegativeBinomial () :
		TypedPrior<negative_binomial_h_param_t, negative_binomial_q_param_t>(),
		a_R(0.0), b_R(0.0) , a_P(0.0), b_P(0.0),
		LSDR_M(0.0) , LSDP_M(0.0) ,
		 a(0.0), b(0.0), lsd(0.0) {};

public:
	PriorNegativeBinomial (double initgamma, double R_M, double P_M , double LSDR_M, double LSDP_M, double a_R, double b_R, double a_P, double b_P , double a, double b, double lsd) :
		TypedPrior<negative_binomial_h_param_t, negative_binomial_q_param_t>(negative_binomial_h_param_t(initgamma), negative_binomial_q_param_t(R_M,P_M),false)
		, a_R(a_R), b_R(b_R) , a_P(a_P), b_P(b_P) ,
		LSDR_M(LSDR_M) , LSDP_M(LSDP_M),
		 a(a), b(a), lsd(lsd)  {
			assert ((P_M <=1) && (P_M >=0));
	}

	void init_q_param () {}
	void init_h_param () {}

	void update_q_param (const  double U, const  int K) ;
	void update_h_param (const double U ,  const int K , const  std::vector<int> & nj );

	int init_M_na() ;
	int update_M_na(const double U ,  const int K);

	double get_gamma() {return this->h_param.gamma;};

};
// --------------------------------------------------------------------------------------------------------------------

#endif /* ANTMAN_SRC_PRIORNEGATIVEBINOMIAL_HPP_ */
