/*
 * GibbsResultRCpp.hpp
 *
 *  Created on: Jun 14, 2019
 *      Author: toky
 */

#ifndef ANTMAN_SRC_GIBBSRESULTRCPP_HPP_
#define ANTMAN_SRC_GIBBSRESULTRCPP_HPP_

#include "GibbsResult.hpp"

class GibbsResultRCpp : public GibbsResult {

private :
	int iteration;

public:
	// Ouput Variables
	std::vector<Rcpp::IntegerVector> CI;
	std::vector<void*>               H;
	std::vector<long>                K;
	std::vector<long>                M;
	std::vector<long>                MNA;
	std::vector<void*>               Q;
	std::vector<arma::vec>           S;
	std::vector<Rcpp::List>          TAU;

	int niter;
	int output_codes;

	GibbsResultRCpp (int niter, int output_codes)  ;

	 void log_output (
			 cluster_indices_t& ci_current,
			 arma::vec & S_current,
			 unsigned int M,
			 unsigned int K,
			 unsigned int M_na,
			 Mixture * mixture,
			 Prior * prior) ;


		Rcpp::List getList () ;
};


#endif /* ANTMAN_SRC_GIBBSRESULTRCPP_HPP_ */
