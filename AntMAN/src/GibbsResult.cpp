/*
 * GibbsResult.cpp
 *
 *  Created on: Jun 13, 2019
 *      Author: toky
 */

#include <RcppArmadillo.h>
#include "GibbsResult.hpp"


#include "Prior.hpp"
#include "Mixture.hpp"


	GibbsResult::GibbsResult(int niter, int output_codes) : iteration(0),  niter(niter), output_codes (output_codes) {

		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_CI )) { VERBOSE_INFO("Record    CI   "); CI.resize(niter);   }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_TAU)) { VERBOSE_INFO("Record    TAU  "); TAU.resize(niter);  }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_S  )) { VERBOSE_INFO("Record    S    "); S.resize(niter);    }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_M  )) { VERBOSE_INFO("Record    M    "); M.resize(niter);    }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_K  )) { VERBOSE_INFO("Record    K    "); K.resize(niter);    }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_Mna)) { VERBOSE_INFO("Record    Mna  "); Mna.resize(niter);  }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_H  )) { VERBOSE_INFO("Record    H    "); H.resize(niter);    }
		if (AM_OUTPUT_HAS(output_codes,AM_OUTPUT_Q  )) { VERBOSE_INFO("Record    Q    "); Q.resize(niter);    }
	}

	 void GibbsResult::log_output (
			 cluster_indices_t& ci_current,
			 Rcpp::NumericVector & S_current,
			 unsigned int M,
			 unsigned int K,
			 unsigned int M_na,
			 Mixture * mixture,
			 Prior * prior) {



		 GibbsResult & result = *this;
		 int output           = this->output_codes;

		 if (AM_OUTPUT_HAS(output,AM_OUTPUT_CI)) {
			 result.CI[iteration]=ci_current;
		 };

		 if (AM_OUTPUT_HAS(output,AM_OUTPUT_S)) {
			 result.S[iteration]=S_current;
		 }

		 if (AM_OUTPUT_HAS(output,AM_OUTPUT_M)) {
			 result.M[iteration]=M;
		 }

		if (AM_OUTPUT_HAS(output,AM_OUTPUT_K))  {
			result.K[iteration]=K;
		}

		if (AM_OUTPUT_HAS(output,AM_OUTPUT_Mna)) {
			result.Mna[iteration]=M_na;
		}

		if (AM_OUTPUT_HAS(output,AM_OUTPUT_TAU)) {
			result.TAU[iteration]=mixture->get_tau();
		}

		if (AM_OUTPUT_HAS(output,AM_OUTPUT_H)) {
			VERBOSE_ERROR("Unsupported case: AM_OUTPUT_H");
		}

		if (AM_OUTPUT_HAS(output,AM_OUTPUT_Q)) {
			VERBOSE_ERROR("Unsupported case: AM_OUTPUT_Q");
		}

		iteration++;
	}


		Rcpp::List GibbsResult::getList () {

			Rcpp::List list = Rcpp::List::create(
					Rcpp::Named("CI")    =  this->CI ,
					Rcpp::Named("S")     =  this->S  ,
					Rcpp::Named("M")     =  this->M  ,
					Rcpp::Named("K")     =  this->K  ,
					Rcpp::Named("Mna")   =  this->Mna,
					Rcpp::Named("Tau")   =  this->TAU
								);


			//					Rcpp::Named("H")     =  res.H  ,
			//					Rcpp::Named("Q")     =  res.Q
			//					Rcpp::Named("TAU")   =  res.TAU,

			 if (not this->CI .size()) list["CI"]  = R_NilValue ;
			 if (not this->TAU.size()) list["Tau"] = R_NilValue ;
			 if (not this->S  .size()) list["S"]   = R_NilValue ;
			 if (not this->M  .size()) list["M"]   = R_NilValue ;
			 if (not this->K  .size()) list["K"]   = R_NilValue ;
			 if (not this->Mna.size()) list["Mna"] = R_NilValue ;
			 if (not this->H  .size()) list["H"]   = R_NilValue ;
			 if (not this->Q  .size()) list["Q"]   = R_NilValue ;

			 return list;


		}



