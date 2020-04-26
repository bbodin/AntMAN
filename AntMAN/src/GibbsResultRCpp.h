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
		{"CI" ,           {"Clusters allocation", 1 << 0}},
		{"TAU" ,          {"", 1 << 1}},
		{"W" ,            {"", 1 << 2}},
		{"U" ,            {"", 1 << 3}},
		{"M" ,            {"", 1 << 4}},
		{"K" ,            {"Number of clusters", 1 << 5}},
		{"H" ,            {"", 1 << 7}},
		{"Q" ,            {"", 1 << 8}},
		{"PREDICTIVE" ,   {"", 1 << 9}},
		{"ALL" ,          {"", (1 << 10) - 1}},
};


static inline bool AM_OUTPUT_HAS (int CODE, int ITEM) {
	return CODE & ITEM;
}

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
	std::vector<double>              U;
	std::vector<long>                M;
	std::vector<Rcpp::List>          Q;
	std::vector<arma::vec>           W;
	std::vector<arma::vec>           PREDICTIVE;
	std::vector<Rcpp::List>          TAU;

	int niter;
	int output_codes;

	GibbsResultRCpp (int niter, int output_codes)  ;

	 void log_output (
			 const cluster_indices_t& CI,
			 const arma::vec & W,
			 const arma::vec & PREDICTIVE,
			 const double U,
			 const unsigned int M,
			 const unsigned int K,
			 const Mixture * mixture,
			 const Prior * prior) ;


		Rcpp::List getList () ;
};


#endif /* ANTMAN_SRC_GIBBSRESULTRCPP_H_ */
