/*
 * GibbsResult.hpp
 *
 *  Created on: Apr 1, 2019
 *      Author: toky
 */

#ifndef ANTMAN_SRC_GIBBSRESULT_HPP_
#define ANTMAN_SRC_GIBBSRESULT_HPP_




struct GibbsResult {

	// Ouput Variables
	std::vector<double> lambda;
	std::vector<double> gamma;
	std::vector<double> U;

	std::vector<Rcpp::IntegerVector> ci;
	std::vector<Rcpp::NumericVector> S;
	std::vector<long> M;
	std::vector<long> K;

};


#endif /* ANTMAN_SRC_GIBBSRESULT_HPP_ */
