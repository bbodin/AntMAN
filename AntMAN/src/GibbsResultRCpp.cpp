/*
 * GibbsResult.cpp
 *
 *  Created on: Jun 13, 2019
 *      Author: toky
 */

#include <RcppArmadillo.h>
#include "GibbsResultRCpp.hpp"


#include "Prior.hpp"
#include "Mixture.hpp"


	GibbsResultRCpp::GibbsResultRCpp(int niter, int output_codes) : iteration(0),  niter(niter), output_codes (output_codes) {

		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("CI").code)) { VERBOSE_INFO("Record    CI   "); CI.resize(niter);   }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("TAU").code)) { VERBOSE_INFO("Record    TAU  "); TAU.resize(niter);  }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("S").code) ){ VERBOSE_INFO("Record    S    "); S.resize(niter);    }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("M").code) ){ VERBOSE_INFO("Record    M    "); M.resize(niter);    }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("K").code) ){ VERBOSE_INFO("Record    K    "); K.resize(niter);    }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("MNA").code)) { VERBOSE_INFO("Record    Mna  "); MNA.resize(niter);  }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("H").code) ){ VERBOSE_INFO("Record    H    "); H.resize(niter);    }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("Q").code) ){ VERBOSE_INFO("Record    Q    "); Q.resize(niter);    }
	}

	 void GibbsResultRCpp::log_output (
			 cluster_indices_t& ci_current,
			 arma::vec & S_current,
			 unsigned int M,
			 unsigned int K,
			 unsigned int M_na,
			 Mixture * mixture,
			 Prior * prior) {

		 GibbsResultRCpp & result = *this;

		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("CI").code)) {result.CI[iteration]=ci_current;};
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("S").code)) {result.S[iteration]=S_current;}
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("M").code) ) {result.M[iteration]=M;		 }
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("K").code) ) {result.K[iteration]=K;		}
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("MNA").code) ) {result.MNA[iteration]=M_na;		}
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("TAU").code)) {result.TAU[iteration]=mixture->get_tau();		}
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("H").code) ) {result.H[iteration]=prior->get_h()->get_Rcpp_list();	}
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("Q").code) ) {result.H[iteration]=prior->get_q()->get_Rcpp_list();	}

		iteration++;
	}


		Rcpp::List GibbsResultRCpp::getList () {

			Rcpp::List list = Rcpp::List::create(
					Rcpp::Named("CI")    =  this->CI ,
					Rcpp::Named("S")     =  this->S  ,
					Rcpp::Named("M")     =  this->M  ,
					Rcpp::Named("K")     =  this->K  ,
					Rcpp::Named("MNA")   =  this->MNA,
					Rcpp::Named("TAU")   =  this->TAU,
					Rcpp::Named("H")     =  this->H  ,
					Rcpp::Named("Q")     =  this->Q );

			 if (not this->CI .size()) list["CI"]  = R_NilValue ;
			 if (not this->TAU.size()) list["TAU"] = R_NilValue ;
			 if (not this->S  .size()) list["S"]   = R_NilValue ;
			 if (not this->M  .size()) list["M"]   = R_NilValue ;
			 if (not this->K  .size()) list["K"]   = R_NilValue ;
			 if (not this->MNA.size()) list["MNA"] = R_NilValue ;
			 if (not this->H  .size()) list["H"]   = R_NilValue ;
			 if (not this->Q  .size()) list["Q"]   = R_NilValue ;

			 return list;


		}



