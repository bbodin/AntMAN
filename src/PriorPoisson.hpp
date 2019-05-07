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

struct poisson_q_param_t {
	double lambda ;
	poisson_q_param_t () : lambda(0.0) {} // TODO remove 0.0
	poisson_q_param_t (double lambda) : lambda (lambda) {}
};

struct poisson_h_param_t {
	double gamma ;
	poisson_h_param_t () : gamma(0.0) {} // TODO remove 0.0
	poisson_h_param_t (double gamma) : gamma (gamma) {}
};

// --------------------------------------------------------------------------------------------------------------------

class PriorPoisson : public TypedPrior < poisson_h_param_t, poisson_q_param_t> {

private:

	double ah,bh; //hyper-prior parameters for h
	double aq,bq; //hyper-prior parameters for q
	double lsd; //this is the standard deviation of the MH algorithm to update gamma.

	PriorPoisson () :
		TypedPrior<poisson_h_param_t, poisson_q_param_t>(0.0, 0.0),
		ah(0.0), bh(0.0) ,
		aq(0.0), bq(0.0), lsd(0.0) {};

public:
	PriorPoisson (double inith, double initq, double ah, double bh , double aq, double bq, double lsd) :
		TypedPrior<poisson_h_param_t, poisson_q_param_t>(inith,initq,false),
		ah(ah), bh(bh) ,
		aq(aq), bq(bq), lsd(lsd)   {
			// TODO : requirements ?
	}

	PriorPoisson (double lambda, double gamma) :
		TypedPrior<poisson_h_param_t, poisson_q_param_t>(lambda, gamma),
		ah(0.0), bh(0.0) ,
		aq(0.0), bq(0.0), lsd(0.0)   {
			// TODO : requirements ?
	}

	double get_gamma() {return this->h_param.gamma;};

	void init_q_param () {}
	void init_h_param () {}

	void update_q_param (const double U, const int K) ;
	void update_h_param (const double U, const int K, const std::vector<int> &nj );

	int init_M_na() ;
	int update_M_na(const double U ,  const int K);


};
// --------------------------------------------------------------------------------------------------------------------




#endif /* ANTMAN_SRC_PRIORPOISSON_HPP_ */
