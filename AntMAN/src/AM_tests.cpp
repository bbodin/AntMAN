/*
 * test-PriorNegativeBinomial.cpp
 *
 *  Created on: 16 May 2019
 *      Author: toky
 */



#include "Priors.hpp"

void test_prior (Prior & to_test) {

	std::cout << "run init_M_na : " << std::endl;
	int Mna0 = to_test.init_M_na(1);

	std::cout << " - init_M_na(K=1) = " << Mna0 << std::endl;


	std::cout << "run update_M_na" << std::endl;
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
		std::cout << " - update_M_na(U= " << U  << ",    K=1) = "         << to_test.update_M_na(U,1) << std::endl;

	}



}

PriorPoissonGamma prepare_poisson () {
	PriorPoissonGamma to_test(poisson_gamma_h_param_t(2,1,1,0),poisson_gamma_q_param_t(1));
	return to_test;
}

PriorNegativeBinomial prepare_negative_binomial () {

	negbin_component R;
	negbin_component P;
	R.M = 1;
	P.M = 0.9;
	R.fixed = true;
	P.fixed = true;

    PriorNegativeBinomial to_test(negative_binomial_gamma_h_param_t(2,1,1,0),negative_binomial_gamma_q_param_t(R,P));
	return to_test;
}

// [[Rcpp::export]]
void test_priors () {

	VERBOSE_LEVEL = DEBUG_LEVEL ;

	PriorNegativeBinomial negbin  = prepare_negative_binomial ();
	PriorPoissonGamma     pois    = prepare_poisson ();

	test_prior (pois);
	test_prior (negbin);



}
