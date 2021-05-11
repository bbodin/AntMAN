/*
 * BasicTest.cpp
 *
 *  Created on: May 11, 2021
 *      Author: toky
 */

#define BOOST_TEST_MODULE BasicTest
#include "testutils.h"

#include <Prior.h>
#include <PriorPoisson.h>
#include <PriorNegativeBinomial.h>

BOOST_AUTO_TEST_SUITE( test_suite_prior )


void test_prior (Prior & to_test) {

	VERBOSE_DEBUG("run init_M_na : ");
	int Mna0 = to_test.init_M_na(1);

	VERBOSE_DEBUG(" - init_M_na(K=1) = " << Mna0);


	VERBOSE_DEBUG("run update_M_na");

	std::vector<double> values = {
			335.185,
			0.0325548,
			0.972901,
			1.69605,
			0.00172949,
			1.00167,
			0.0556233,
			0.605071,
			3.35195e-06,
			0.813054,

	};

	for (auto U : values) {
		VERBOSE_DEBUG(" - update_M_na(U= " << U  << ",    K=1) = "         << to_test.update_M_na(U,1));
	}

}


BOOST_AUTO_TEST_CASE( test_poisson_prior )
{
	PriorPoisson to_test(poisson_gamma_h_param_t(2,1,1,0),poisson_gamma_q_param_t(1));
	test_prior (to_test);
}


BOOST_AUTO_TEST_CASE( test_negbin_prior )
{
	negbin_component R;
	negbin_component P;
	R.value = 1;
	P.value = 0.9;
	R.fixed = true;
	P.fixed = true;

	PriorNegativeBinomial to_test(negative_binomial_gamma_h_param_t(2,1,1,0),negative_binomial_gamma_q_param_t(R,P));

	test_prior (to_test);
}

BOOST_AUTO_TEST_SUITE_END()



