/*
 *  AntMAN Package
 *  Unit test
 *
 */



#include "testutils.h"

#include <utils.h>


std::string getString (AntMANLogger& logger) {
	 std::ostringstream outputstr;
	std::vector<std::string> names = {};
	int total = 0;

	if (logger.haslog("CI")) {total++;names.push_back("CI");}
	if (logger.haslog("CI")) {
		outputstr << "Size of CI vector = " << logger.getlog<cluster_indices_t>("CI").size() << std::endl ;
		outputstr << "Last CI vector = " << logger.getlog<cluster_indices_t>("CI").back() << std::endl ;
	}
	auto K_values = logger.getlog<int>("K");

	if (logger.haslog("K")) {
		outputstr << "K = " ;
		for (arma::uword i=0; i < K_values.size(); ++i) {
			outputstr << K_values[i] << ' ';
		}
		outputstr <<  std::endl ;
	}


	return outputstr.str();



}



