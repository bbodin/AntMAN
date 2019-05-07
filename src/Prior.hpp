/*
 * Prior.hpp
 *
 *  Created on: Apr 8, 2019
 *      Author: toky
 */

#ifndef ANTMAN_SRC_PRIOR_HPP_
#define ANTMAN_SRC_PRIOR_HPP_

#include <cassert>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// --------------------------------------------------------------------------------------------------------------------

class Prior {


public :
	virtual void    update (const double U, const int K, const std::vector<int> &nj ) = 0 ;

	virtual double  get_gamma() const = 0;

	virtual int     init_M_na()= 0;
	virtual int     update_M_na(const double U ,  const int K)= 0;

	virtual        ~Prior() {};

};

template <typename h_param_t, typename q_param_t>
class TypedPrior : public Prior {

protected :
	h_param_t h_param;
	q_param_t q_param;
public:

	void update (const double U, const int K, const std::vector<int> &nj ) {
		 this->q_param.update (U, K, this->h_param);
		 this->h_param.update (U, K, nj, this->q_param);
	};

	double get_gamma() const {return this->h_param.gamma;};

	TypedPrior(h_param_t h_param, q_param_t q_param)             : h_param(h_param), q_param(q_param) {};


	virtual    ~TypedPrior() {};

};


#endif /* ANTMAN_SRC_PRIOR_HPP_ */
