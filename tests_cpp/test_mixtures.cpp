/*
 *  AntMAN Package
 *
 */

#ifndef NO_RCPP
#define NO_RCPP
#endif

#include <iostream>
#include <boost/program_options.hpp>
#include "../AntMAN/src/verbose.h"
#include "AM_test_mixtures.h"
#include "AM_test_priors.h"

namespace po = boost::program_options;

int main (int ac, char** av) {

	po::options_description desc("Allowed options");
	desc.add_options()
	    ("help", "Help message")
	    ("verbose", po::value<int>()->default_value(LOG_LEVEL), "Verbosity level")
	    ("niter", po::value<long>()->default_value(5000), "Total number of iteration")
	    ("burnin", po::value<long>()->default_value(1000), "Number of iteration to drop")
	    ("thin", po::value<long>()->default_value(10), "Rate of storing after burnin")
	    ("input-file", po::value< std::vector<std::string> >(), "input file")
	;



	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);

	if (vm.count("verbose")) {
		VERBOSE_LEVEL(vm["verbose"].as<int>());
	}

	VERBOSE_LOG("VERBOSE_LEVEL=" << VERBOSE_LEVEL());
	VERBOSE_INFO("INFO_LEVEL ACTIVATED");
	VERBOSE_DEBUG("DEBUG_LEVEL ACTIVATED");

	if (vm.count("help")) {
		COUT_STREAM << desc << "\n";
	    return 0;
	}
	long thin =  vm["thin"].as<long>();
	long niter =  vm["niter"].as<long>();
	long burnin =  vm["burnin"].as<long>();

	test_mixtures (niter ,burnin ,thin);
}