/*
 * GibbsResultsToFile.cpp
 *
 *  Created on: Aug 2, 2019
 */


#include "GibbsResult.h"

void GibbsResultIntoFile ::log_output (
			 cluster_indices_t& ci_current,
			 arma::vec & S_current,
			 unsigned int M,
			 unsigned int K,
			 unsigned int M_na,
			 Mixture * mixture,
			 Prior * prior) {

	std::string header = this->_dirname + "/" + "AntMan_";





	//write each file
	std::ofstream fdesc;
	fdesc.open (header + "CI.txt");	for (auto idx : ci_current) {fdesc << idx;} fdesc << std::endl; fdesc.close();
	fdesc.open (header + "M.txt");	fdesc << M; fdesc << std::endl; fdesc.close();
	fdesc.open (header + "K.txt");	fdesc << K; fdesc << std::endl; fdesc.close();
	fdesc.open (header + "Mna.txt");	fdesc << M_na; fdesc << std::endl; fdesc.close();
	fdesc.open (header + "S.txt");	for (auto idx : S_current) {fdesc << idx;} fdesc << std::endl; fdesc.close();
	//fdesc.open (header + "TAU.txt");	fdesc << mixture->get_tau(); fdesc << std::endl; fdesc.close();
	//fdesc.open (header + "H.txt");	fdesc << prior->get_h()->get_str(); fdesc << std::endl; fdesc.close();
	//fdesc.open (header + "Q.txt");	fdesc << prior->get_q()->get_str(); fdesc << std::endl; fdesc.close();





}
