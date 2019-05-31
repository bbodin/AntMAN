/*
 * test-PriorNegativeBinomial.cpp
 *
 *  Created on: 16 May 2019
 *      Author: toky
 */



#include <PriorNegativeBinomial.hpp>

int main () {

	std::cout << "test starts" << std::endl;

	std::cout << "init components" << std::endl;


	negbin_component R;
	negbin_component P;
	R.M = 1;
	P.M = 0.1;
	R.fixed = true;
	P.fixed = true;


	std::cout << "init H and Q" << std::endl;


	negative_binomial_gamma_q_param_t * q = new negative_binomial_gamma_q_param_t(R,P);
	negative_binomial_gamma_h_param_t * h = new negative_binomial_gamma_h_param_t(2,1,1,0);

	std::cout << "init PNB" << std::endl;


	PriorNegativeBinomial to_test(*h,*q);
	std::cout << "run init Mna" << std::endl;


	int Mna0 = to_test.init_M_na(1000);

	VERBOSE_INFO("Mna0 = " << Mna0);

	return 0 ;

}
