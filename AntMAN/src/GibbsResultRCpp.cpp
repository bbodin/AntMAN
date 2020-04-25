/*
 * GibbsResult.cpp
 *
 *  Created on: Jun 13, 2019
 */

#include "GibbsResultRCpp.h"
#include "Mixture.h"
#include "Prior.h"


	GibbsResultRCpp::GibbsResultRCpp(int niter, int output_codes) : iteration(0),  niter(niter), output_codes (output_codes) {

		for (auto idx : AM_OUTPUTS){
			if (AM_OUTPUT_HAS(output_codes,idx.second.code)) {
				VERBOSE_INFO("Record " << idx.second.name);
			}
		}
	}

	 void GibbsResultRCpp::log_output (
			 cluster_indices_t& CI,
			 arma::vec&         W,
			 arma::vec&         PREDICTIVE,
			 double             U,
			 unsigned int       M,
			 unsigned int       K,
			 Mixture *          mixture,
			 Prior *            prior) {

		 GibbsResultRCpp & result = *this;

		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("CI").code)) {result.CI.push_back(CI);}
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("W").code)) {result.W.push_back(W);}
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("PREDICTIVE").code)) {result.PREDICTIVE.push_back(PREDICTIVE);}
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("U").code) ) {result.U.push_back(U);		 }
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("M").code) ) {result.M.push_back(M);		 }
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("K").code) ) {result.K.push_back(K);		}
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("TAU").code)) {result.TAU.push_back(mixture->get_tau());		}
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("H").code) ) {result.H.push_back(prior->get_h()->get_Rcpp_list());	}
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("Q").code) ) {result.Q.push_back(prior->get_q()->get_Rcpp_list());	}

		iteration++;
	}


		Rcpp::List GibbsResultRCpp::getList () {

			std::vector<std::string> names = {};
			int total = 0;
			if (this->CI .size()) {total++;names.push_back("CI");}
			if (this->TAU.size()) {total++;names.push_back("TAU");}
			if (this->W  .size()) {total++;names.push_back("W");}
			if (this->PREDICTIVE  .size()) {total++;names.push_back("PREDICTIVE");}
			if (this->U  .size()) {total++;names.push_back("U");}
			if (this->M  .size()) {total++;names.push_back("M");}
			if (this->K  .size()) {total++;names.push_back("K");}
			if (this->H  .size()) {total++;names.push_back("H");}
			if (this->Q  .size()) {total++;names.push_back("Q");}

			Rcpp::List my_list(total);
			my_list.attr("names") = names;

			int cnt = 0;
			if (this->CI .size()) {my_list[cnt++] =  this->CI ;}
			if (this->TAU.size()) {my_list[cnt++] =  this->TAU;}
			if (this->W  .size()) {my_list[cnt++] =  this->W  ;}
			if (this->PREDICTIVE  .size()) {my_list[cnt++] =  this->PREDICTIVE  ;}
			if (this->U  .size()) {my_list[cnt++] =  this->U  ;}
			if (this->M  .size()) {my_list[cnt++] =  this->M  ;}
			if (this->K  .size()) {my_list[cnt++] =  this->K  ;}
			if (this->H  .size()) {my_list[cnt++] =  this->H  ;}
			if (this->Q  .size()) {my_list[cnt++] =  this->Q  ;}

			return my_list;



		}



