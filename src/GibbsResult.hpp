/*
 * GibbsResult.hpp
 *
 *  Created on: Apr 1, 2019
 *      Author: toky
 */

#ifndef ANTMAN_SRC_GIBBSRESULT_HPP_
#define ANTMAN_SRC_GIBBSRESULT_HPP_


static const std::string  AM_OUTPUT_STR_CI  = "CI";
static const std::string  AM_OUTPUT_STR_TAU = "Tau";
static const std::string  AM_OUTPUT_STR_S   = "S";
static const std::string  AM_OUTPUT_STR_M   = "M";
static const std::string  AM_OUTPUT_STR_K   = "K";
static const std::string  AM_OUTPUT_STR_Mna = "Mna";
static const std::string  AM_OUTPUT_STR_H   = "H";
static const std::string  AM_OUTPUT_STR_Q   = "Q";
static const std::string  AM_OUTPUT_STR_ALL = "all";


static const int  AM_OUTPUT_CI  = 1 << 0;
static const int  AM_OUTPUT_TAU = 1 << 1;
static const int  AM_OUTPUT_S   = 1 << 2;
static const int  AM_OUTPUT_M   = 1 << 3;
static const int  AM_OUTPUT_K   = 1 << 4;
static const int  AM_OUTPUT_Mna = 1 << 5;
static const int  AM_OUTPUT_H   = 1 << 6;
static const int  AM_OUTPUT_Q   = 1 << 7;

static const int  AM_OUTPUT_DEFAULT   = AM_OUTPUT_CI | AM_OUTPUT_S | AM_OUTPUT_M | AM_OUTPUT_K;
static const int  AM_OUTPUT_ALL       = AM_OUTPUT_CI | AM_OUTPUT_TAU | AM_OUTPUT_S | AM_OUTPUT_M | AM_OUTPUT_K | AM_OUTPUT_Mna | AM_OUTPUT_H | AM_OUTPUT_Q;

static inline bool AM_OUTPUT_HAS (int CODE, int ITEM) {return CODE & ITEM;};

static inline int AM_GENERATOR_OUTPUT_CODE (std::vector <std::string> output) {
		int res = 0 ;
		for (std::string e : output ) {
			     if (e == AM_OUTPUT_STR_ALL) { VERBOSE_DEBUG("Add all "); return  AM_OUTPUT_ALL ; }
			else if (e == AM_OUTPUT_STR_CI ) { VERBOSE_DEBUG("Add CI  "); res |=  AM_OUTPUT_CI  ; }
			else if (e == AM_OUTPUT_STR_TAU) { VERBOSE_DEBUG("Add TAU "); res |=  AM_OUTPUT_TAU ; }
			else if (e == AM_OUTPUT_STR_S  ) { VERBOSE_DEBUG("Add S   "); res |=  AM_OUTPUT_S   ; }
			else if (e == AM_OUTPUT_STR_M  ) { VERBOSE_DEBUG("Add M   "); res |=  AM_OUTPUT_M   ; }
			else if (e == AM_OUTPUT_STR_K  ) { VERBOSE_DEBUG("Add K   "); res |=  AM_OUTPUT_K   ; }
			else if (e == AM_OUTPUT_STR_Mna) { VERBOSE_DEBUG("Add Mna "); res |=  AM_OUTPUT_Mna ; }
			else if (e == AM_OUTPUT_STR_H  ) { VERBOSE_DEBUG("Add H   "); res |=  AM_OUTPUT_H   ; }
			else if (e == AM_OUTPUT_STR_Q  ) { VERBOSE_DEBUG("Add Q   "); res |=  AM_OUTPUT_Q   ; }
			else {VERBOSE_ERROR ("Unsupported output item: " << e);}
		}
		return res;
	};



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

	int niter;
	int output_codes;

	GibbsResult(int niter, int output_codes) : niter(niter), output_codes (output_codes) {

		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_CI )) { VERBOSE_INFO("Record    CI   "); CI.resize(niter);   }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_TAU)) { VERBOSE_INFO("Record    TAU  "); TAU.resize(niter);  }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_S  )) { VERBOSE_INFO("Record    S    "); S.resize(niter);    }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_M  )) { VERBOSE_INFO("Record    M    "); M.resize(niter);    }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_K  )) { VERBOSE_INFO("Record    K    "); K.resize(niter);    }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_Mna)) { VERBOSE_INFO("Record    Mna  "); Mna.resize(niter);  }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_H  )) { VERBOSE_INFO("Record    H    "); H.resize(niter);    }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_Q  )) { VERBOSE_INFO("Record    Q    "); Q.resize(niter);    }
	}

};


#endif /* ANTMAN_SRC_GIBBSRESULT_HPP_ */
