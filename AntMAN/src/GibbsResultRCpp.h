/*
 * GibbsResultRCpp.hpp
 *
 *  Created on: Jun 14, 2019
 */

#ifndef ANTMAN_SRC_GIBBSRESULTRCPP_H_
#define ANTMAN_SRC_GIBBSRESULTRCPP_H_

#ifndef HAS_RCPP
#assert HAS_RCPP
#undef  NO_RCPP
#define HAS_RCPP
#endif


#include "GibbsResult.h"


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



class GibbsResultRCpp : public GibbsResult {

private :
	int iteration;

public:
	// Ouput Variables
	std::vector<arma::ivec>          CI;
	std::vector<Rcpp::List>          H;
	std::vector<long>                K;
	std::vector<long>                M;
	std::vector<long>                MNA;
	std::vector<Rcpp::List>          Q;
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


#endif /* ANTMAN_SRC_GIBBSRESULTRCPP_H_ */
