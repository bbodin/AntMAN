/*
 * GibbsResult.hpp
 *
 *  Created on: Apr 1, 2019
 *      Author: toky
 */

#ifndef ANTMAN_SRC_GIBBSRESULT_HPP_
#define ANTMAN_SRC_GIBBSRESULT_HPP_

#include "utils.hpp"

class Mixture;
class Prior;

struct Output_Type {
	std::string name;
	int         code;
};

static const std::map<std::string, Output_Type> AM_OUTPUTS = {
		{"CI" ,  {"Clusters allocation", 1 << 0}},
		{"TAU" , {"", 1 << 1}},
		{"S" ,   {"", 1 << 2}},
		{"M" ,   {"", 1 << 3}},
		{"K" ,   {"Number of clusters", 1 << 4}},
		{"MNA" , {"", 1 << 5}},
		{"H" ,   {"", 1 << 6}},
		{"Q" ,   {"", 1 << 7}},
		{"ALL" , {"", (1 << 8) - 1}},
};


static inline bool AM_OUTPUT_HAS (int CODE, int ITEM) {return CODE & ITEM;}

static inline int AM_GENERATOR_OUTPUT_CODE (std::vector <std::string> output) {
		int res = 0 ;
		for (std::string e : output ) {
			transform(e.begin(), e.end(), e.begin(), ::toupper);
			bool unchanged = true;
			for (std::pair<std::string, Output_Type> map_idx : AM_OUTPUTS) {
			    	 if (e == map_idx.first ) {
			    		 VERBOSE_DEBUG("Add " <<map_idx.second.name); res |=  map_idx.second.code  ;
			    		 unchanged = false;
			    	 }
			}

			if (unchanged) VERBOSE_WARNING ("Unsupported output item: " << e);
		}
		return res;
	}


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
public:
	GibbsResultPlain () {} ;

	 void log_output (
			 cluster_indices_t& ci_current,
			 arma::vec & S_current,
			 unsigned int M,
			 unsigned int K,
			 unsigned int M_na,
			 Mixture * mixture,
			 Prior * prior) {} ;

};


#endif /* ANTMAN_SRC_GIBBSRESULT_HPP_ */
