/*
 * GibbsResult.hpp
 *
 *  Created on: Apr 1, 2019
 */

#ifndef ANTMAN_SRC_GIBBSRESULT_H_
#define ANTMAN_SRC_GIBBSRESULT_H_

#include "utils.h"

class Mixture;
class Prior;

class GibbsResult {
public :
	virtual ~GibbsResult() {} ;

	virtual void log_output (
				 cluster_indices_t& ci_current,
				 arma::vec & S_current,
				 unsigned int M,
				 unsigned int K,
				 unsigned int M_na,
				 Mixture * mixture,
				 Prior * prior) = 0;
};


class GibbsResultPlain : public GibbsResult {
private :
	long int index;
	arma::vec Ks;
public:
	arma::vec& K () {return Ks;};
	GibbsResultPlain (long int maxsize) : index (0) ,  Ks (maxsize){} ;

	 void log_output (
			 cluster_indices_t& ci_current,
			 arma::vec & S_current,
			 unsigned int M,
			 unsigned int K,
			 unsigned int M_na,
			 Mixture * mixture,
			 Prior * prior) {
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
			 cluster_indices_t& ci_current,
			 arma::vec & S_current,
			 unsigned int M,
			 unsigned int K,
			 unsigned int M_na,
			 Mixture * mixture,
			 Prior * prior) ;

};


#endif /* ANTMAN_SRC_GIBBSRESULT_H_ */
