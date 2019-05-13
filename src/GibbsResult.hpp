/*
 * GibbsResult.hpp
 *
 *  Created on: Apr 1, 2019
 *      Author: toky
 */

#ifndef ANTMAN_SRC_GIBBSRESULT_HPP_
#define ANTMAN_SRC_GIBBSRESULT_HPP_


static const int  AM_OUTPUT_CI  = 1 << 0;
static const int  AM_OUTPUT_TAU = 1 << 1;
static const int  AM_OUTPUT_S   = 1 << 2;
static const int  AM_OUTPUT_M   = 1 << 3;
static const int  AM_OUTPUT_K   = 1 << 4;
static const int  AM_OUTPUT_Mna = 1 << 5;
static const int  AM_OUTPUT_H   = 1 << 6;
static const int  AM_OUTPUT_Q   = 1 << 7;

static const int  AM_OUTPUT_DEFAULT   = AM_OUTPUT_CI & AM_OUTPUT_S & AM_OUTPUT_M & AM_OUTPUT_K;

static inline bool AM_OUTPUT_HAS (int CODE, int ITEM) {return CODE && ITEM;};


class GibbsResult {
public:
	// Ouput Variables
	std::vector<Rcpp::IntegerVector> CI;
	std::vector<void*>               TAU;
	std::vector<Rcpp::NumericVector> S;
	std::vector<long>                M;
	std::vector<long>                K;
	std::vector<long>                Mna;
	std::vector<void*>               H;
	std::vector<void*>               Q;

	GibbsResult(int niter) {
		CI.resize(niter);
		TAU.resize(niter);
		S.resize(niter);
		M.resize(niter);
		K.resize(niter);
		Mna.resize(niter);
		H.resize(niter);
		Q.resize(niter);
	}

};


#endif /* ANTMAN_SRC_GIBBSRESULT_HPP_ */
