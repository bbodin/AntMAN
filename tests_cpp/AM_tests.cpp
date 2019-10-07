/*
 * AM_tests.cpp
 *
 *  Created on: Jun 14, 2019
 *      Author: toky
 */

#include <iostream>
#include <boost/program_options.hpp>
#include "../AntMAN/src/verbose.h"

namespace po = boost::program_options;

void test_priors   ();
void test_mixtures (long,long,long);

int main (int ac, char** av) {

	po::options_description desc("Allowed options");
	desc.add_options()
	    ("help", "Help message")
	    ("verbose", po::value<int>(&VERBOSE_LEVEL)->default_value(LOG_LEVEL), "Verbosity level")
	    ("niter", po::value<long>()->default_value(5000), "Total number of iteration")
	    ("burnin", po::value<long>()->default_value(1000), "Number of iteration to drop")
	    ("thin", po::value<long>()->default_value(10), "Rate of storing after burnin")
	    ("input-file", po::value< std::vector<std::string> >(), "input file")
	;



	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);
	VERBOSE_LOG("VERBOSE_LEVEL=" << VERBOSE_LEVEL);

	if (vm.count("help")) {
		COUT_STREAM << desc << "\n";
	    return 0;
	}
	long thin =  vm["thin"].as<long>();
	long niter =  vm["niter"].as<long>();
	long burnin =  vm["burnin"].as<long>();

	test_mixtures (niter ,burnin ,thin);
}
