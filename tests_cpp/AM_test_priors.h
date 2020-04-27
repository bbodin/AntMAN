/*
 * AM_test_priors.h
 *
 *  Created on: Mar 6, 2020
 *      Author: toky
 */


#ifndef TESTS_CPP_AM_TEST_PRIORS_H_
#define TESTS_CPP_AM_TEST_PRIORS_H_

#ifndef NO_RCPP
#define NO_RCPP
#endif

#include "../AntMAN/src/Priors.h"

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

PriorPoisson prepare_poisson () {
	PriorPoisson to_test(poisson_gamma_h_param_t(2,1,1,0),poisson_gamma_q_param_t(1));
	return to_test;
}

PriorNegativeBinomial prepare_negative_binomial () {

	negbin_component R;
	negbin_component P;
	R.value = 1;
	P.value = 0.9;
	R.fixed = true;
	P.fixed = true;

    PriorNegativeBinomial to_test(negative_binomial_gamma_h_param_t(2,1,1,0),negative_binomial_gamma_q_param_t(R,P));
	return to_test;
}

void test_priors () {

	VERBOSE_LEVEL ( LOG_LEVEL );

	PriorNegativeBinomial negbin  = prepare_negative_binomial ();
	PriorPoisson     pois    = prepare_poisson ();

	test_prior (pois);
	test_prior (negbin);



}




#endif /* TESTS_CPP_AM_TEST_PRIORS_H_ */
