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
			if (AM_OUTPUT_HAS(output_codes,idx.second.code)) { VERBOSE_INFO("Record " << idx.second.name); }
		}
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

		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("CI").code)) {result.CI.push_back(ci_current);}
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("S").code)) {result.S.push_back(S_current);}
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("M").code) ) {result.M.push_back(M);		 }
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("K").code) ) {result.K.push_back(K);		}
		 if (AM_OUTPUT_HAS(output_codes,AM_OUTPUTS.at("MNA").code) ) {result.MNA.push_back(M_na);		}
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
			if (this->S  .size()) {total++;names.push_back("S");}
			if (this->M  .size()) {total++;names.push_back("M");}
			if (this->K  .size()) {total++;names.push_back("K");}
			if (this->MNA.size()) {total++;names.push_back("MNA");}
			if (this->H  .size()) {total++;names.push_back("H");}
			if (this->Q  .size()) {total++;names.push_back("Q");}

			Rcpp::List my_list(total);
			my_list.attr("names") = names;

			int cnt = 0;
			if (this->CI .size()) {my_list[cnt++] =  this->CI ;}
			if (this->TAU.size()) {my_list[cnt++] =  this->TAU;}
			if (this->S  .size()) {my_list[cnt++] =  this->S  ;}
			if (this->M  .size()) {my_list[cnt++] =  this->M  ;}
			if (this->K  .size()) {my_list[cnt++] =  this->K  ;}
			if (this->MNA.size()) {my_list[cnt++] =  this->MNA;}
			if (this->H  .size()) {my_list[cnt++] =  this->H  ;}
			if (this->Q  .size()) {my_list[cnt++] =  this->Q  ;}

			return my_list;



		}



