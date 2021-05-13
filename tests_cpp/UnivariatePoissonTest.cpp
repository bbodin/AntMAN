/*
 * UnivariateNormalTest.cpp
 *
 *  Created on: May 11, 2021
 *      Author: toky
 */

#define BOOST_TEST_MODULE UnivariatePoissonTest

#include <boost/test/included/unit_test.hpp>

#include <PriorPoisson.h>
#include <MixtureUnivariatePoisson.h>

void test_Mixture_UnivariatePoisson(long niter, long burnin, long thin) {


	static const arma::ivec y_uvn =  { 9,  9,  9,  9,  9, 10, 10, 16, 16, 18,18,
			18, 18, 19, 19, 19, 19, 19, 19, 19,19, 19, 19, 19,
			19, 19, 19, 19, 19, 19,19, 20, 20, 20, 20, 20, 20,
			20, 20, 20,20, 20, 20, 20, 21, 21, 21, 21, 21, 21,
			22, 22, 22, 22, 22, 22, 22, 22, 22, 22,22, 23, 23,
			23, 23, 23, 23, 23, 23, 23,24, 24, 24, 24, 24, 24, 25, 26,
			26, 32,32, 34 };
	// TODO : Need to go over variable and see what should be checked
	PriorPoisson *prior = new PriorPoisson(poisson_gamma_h_param_t(2,1,1,0.00001),poisson_gamma_q_param_t(3,1,1));
	MixtureUnivariatePoisson * mixture = new MixtureUnivariatePoisson (2, 0.2);
	cluster_indices_t initial_clustering (y_uvn.size());
	AntMANLogger * logger = new AntMANLogger(std::vector<std::string>(), (niter - burnin) / thin );

	auto start_gibbs           = std::chrono::system_clock::now();
	mixture->fit(y_uvn , initial_clustering, false, prior , niter ,burnin ,thin , false , logger);
	auto end_gibbs             = std::chrono::system_clock::now();
	auto elapsed_gibbs         = end_gibbs - start_gibbs;
	auto total_gibbs           = elapsed_gibbs.count() / 1000000.0;
	COUT_STREAM << "Total time: " << total_gibbs << "ms"  << std::endl ;


}


BOOST_AUTO_TEST_SUITE( test_suite_univariate_poisson_dataset )

BOOST_AUTO_TEST_CASE( test_univariate_poisson_100_10_10 )
{
	test_Mixture_UnivariatePoisson(100, 10, 10);
}


BOOST_AUTO_TEST_SUITE_END()




