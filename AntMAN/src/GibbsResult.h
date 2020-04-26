/*
 * GibbsResult.hpp
 *
 *  Created on: Apr 1, 2019
 */

#ifndef ANTMAN_SRC_GIBBSRESULT_H_
#define ANTMAN_SRC_GIBBSRESULT_H_

#include "utils.h"
#include <string>

class Mixture;
class Prior;

class GibbsResult {
public :
	virtual ~GibbsResult() {} ;

	virtual void log_output (
				 const cluster_indices_t& CI,
				 const arma::vec&         W,
				 const arma::vec&         PREDICTIVE,
				 const double             U,
				 const unsigned int       M,
				 const unsigned int       K,
				 const Mixture *          mixture,
				 const Prior *            prior) = 0;


};


class GibbsResultPlain : public GibbsResult {
private :
	long int index;
	arma::vec Ks;
public:
	arma::vec& K () {return Ks;};
	GibbsResultPlain (long int maxsize) : index (0) ,  Ks (maxsize){} ;

	 void log_output (
			 const cluster_indices_t& CI,
			 const arma::vec & W,
			 const arma::vec & PREDICTIVE,
			 const double U,
			 const unsigned int M,
			 const unsigned int K,
			 const Mixture * mixture,
			 const Prior * prior) {
		 Ks[index++] = K;
	 } ;

};

class GibbsResultIntoFile : public GibbsResult {
private :
	size_t count = 0;
	std::string _dirname;
public:
	GibbsResultIntoFile (std::string dirname) : _dirname(dirname) {} ;

	 void log_output (
			 const cluster_indices_t& ci_current,
			 const arma::vec & W,
			 const arma::vec & PREDICTIVE,
			 const double U,
			 const unsigned int M,
			 const unsigned int K,
			 const Mixture * mixture,
			 const Prior * prior) ;

};


#endif /* ANTMAN_SRC_GIBBSRESULT_H_ */
