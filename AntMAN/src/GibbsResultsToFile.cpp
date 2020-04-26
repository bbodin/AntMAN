/*
 * GibbsResultsToFile.cpp
 *
 *  Created on: Aug 2, 2019
 */


#include "GibbsResult.h"

void GibbsResultIntoFile ::log_output (
			 const cluster_indices_t& CI,
			 const arma::vec & W,
			 const arma::vec & PREDICTIVE,
			 const double U,
			 const unsigned int M,
			 const unsigned int K,
			 const Mixture * mixture,
			 const Prior * prior) {

	std::string header = this->_dirname + "/" + "AntMan_";





	//write each file
	std::ofstream fdesc;
	fdesc.open (header + "CI.txt");	for (auto idx : CI) {fdesc << idx;} fdesc << std::endl; fdesc.close();
	fdesc.open (header + "U.txt");	fdesc << U; fdesc << std::endl; fdesc.close();
	fdesc.open (header + "M.txt");	fdesc << M; fdesc << std::endl; fdesc.close();
	fdesc.open (header + "K.txt");	fdesc << K; fdesc << std::endl; fdesc.close();
	fdesc.open (header + "W.txt");	for (auto idx : W) {fdesc << idx;} fdesc << std::endl; fdesc.close();
	//fdesc.open (header + "TAU.txt");	fdesc << mixture->get_tau(); fdesc << std::endl; fdesc.close();
	//fdesc.open (header + "H.txt");	fdesc << prior->get_h()->get_str(); fdesc << std::endl; fdesc.close();
	//fdesc.open (header + "Q.txt");	fdesc << prior->get_q()->get_str(); fdesc << std::endl; fdesc.close();





}
