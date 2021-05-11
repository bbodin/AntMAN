/*
 * MultivariateBinomialTest.cpp
 *
 *  Created on: May 11, 2021
 *      Author: toky
 */

#define BOOST_TEST_MODULE MultivariateBinomialTest
#include "testutils.h"

#include <PriorPoisson.h>
#include <MixtureMultivariateBinomial.h>

void test_Mixture_MultivariateBinomial (long niter, long burnin, long thin) {
	static const arma::imat carcinoma = {
	     //A  B  C  D  E  F  G
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 0, 0, 0 } ,
		 { 0, 0, 0, 0, 1, 0, 0 } ,
		 { 0, 0, 0, 0, 1, 0, 0 } ,
		 { 0, 1, 0, 0, 0, 0, 0 } ,
		 { 0, 1, 0, 0, 0, 0, 0 } ,
		 { 0, 1, 0, 0, 0, 0, 0 } ,
		 { 0, 1, 0, 0, 0, 0, 0 } ,
		 { 0, 1, 0, 0, 0, 0, 0 } ,
		 { 0, 1, 0, 0, 0, 0, 0 } ,
		 { 0, 1, 0, 0, 0, 0, 1 } ,
		 { 0, 1, 0, 0, 1, 0, 0 } ,
		 { 0, 1, 0, 0, 1, 0, 0 } ,
		 { 0, 1, 0, 0, 1, 0, 0 } ,
		 { 0, 1, 0, 0, 1, 0, 0 } ,
		 { 0, 1, 0, 0, 1, 0, 1 } ,
		 { 0, 1, 0, 0, 1, 0, 1 } ,
		 { 0, 1, 0, 0, 1, 0, 1 } ,
		 { 0, 1, 0, 0, 1, 0, 1 } ,
		 { 0, 1, 0, 0, 1, 0, 1 } ,
		 { 1, 0, 0, 0, 0, 0, 0 } ,
		 { 1, 0, 0, 0, 0, 0, 0 } ,
		 { 1, 0, 1, 0, 1, 0, 1 } ,
		 { 1, 1, 0, 0, 0, 0, 0 } ,
		 { 1, 1, 0, 0, 0, 0, 0 } ,
		 { 1, 1, 0, 0, 0, 0, 1 } ,
		 { 1, 1, 0, 0, 1, 0, 0 } ,
		 { 1, 1, 0, 0, 1, 0, 0 } ,
		 { 1, 1, 0, 0, 1, 0, 1 } ,
		 { 1, 1, 0, 0, 1, 0, 1 } ,
		 { 1, 1, 0, 0, 1, 0, 1 } ,
		 { 1, 1, 0, 0, 1, 0, 1 } ,
		 { 1, 1, 0, 0, 1, 0, 1 } ,
		 { 1, 1, 0, 0, 1, 0, 1 } ,
		 { 1, 1, 0, 0, 1, 0, 1 } ,
		 { 1, 1, 0, 0, 1, 1, 1 } ,
		 { 1, 1, 0, 1, 0, 0, 1 } ,
		 { 1, 1, 0, 1, 1, 0, 1 } ,
		 { 1, 1, 0, 1, 1, 0, 1 } ,
		 { 1, 1, 0, 1, 1, 1, 1 } ,
		 { 1, 1, 0, 1, 1, 1, 1 } ,
		 { 1, 1, 0, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 0, 1, 0, 1 } ,
		 { 1, 1, 1, 0, 1, 0, 1 } ,
		 { 1, 1, 1, 0, 1, 0, 1 } ,
		 { 1, 1, 1, 0, 1, 0, 1 } ,
		 { 1, 1, 1, 0, 1, 0, 1 } ,
		 { 1, 1, 1, 0, 1, 0, 1 } ,
		 { 1, 1, 1, 0, 1, 0, 1 } ,
		 { 1, 1, 1, 0, 1, 0, 1 } ,
		 { 1, 1, 1, 0, 1, 0, 1 } ,
		 { 1, 1, 1, 0, 1, 0, 1 } ,
		 { 1, 1, 1, 0, 1, 0, 1 } ,
		 { 1, 1, 1, 0, 1, 0, 1 } ,
		 { 1, 1, 1, 0, 1, 0, 1 } ,
		 { 1, 1, 1, 0, 1, 1, 1 } ,
		 { 1, 1, 1, 0, 1, 1, 1 } ,
		 { 1, 1, 1, 0, 1, 1, 1 } ,
		 { 1, 1, 1, 0, 1, 1, 1 } ,
		 { 1, 1, 1, 0, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 0, 1 } ,
		 { 1, 1, 1, 1, 1, 0, 1 } ,
		 { 1, 1, 1, 1, 1, 0, 1 } ,
		 { 1, 1, 1, 1, 1, 0, 1 } ,
		 { 1, 1, 1, 1, 1, 0, 1 } ,
		 { 1, 1, 1, 1, 1, 0, 1 } ,
		 { 1, 1, 1, 1, 1, 0, 1 } ,
		 { 1, 1, 1, 1, 1, 0, 1 } ,
		 { 1, 1, 1, 1, 1, 0, 1 } ,
		 { 1, 1, 1, 1, 1, 0, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 } ,
		 { 1, 1, 1, 1, 1, 1, 1 }
} ;

int d = carcinoma.n_cols;

arma::vec a0(d);
arma::vec b0(d);
a0.fill(1);
b0.fill(1);

PriorPoisson *priormvb = new PriorPoisson(poisson_gamma_h_param_t(2,1,1,0.00001),poisson_gamma_q_param_t(5,10,2));
MixtureMultivariateBinomial * mixturemvb = new MixtureMultivariateBinomial (a0,b0);
AntMANLogger * logger = new AntMANLogger(std::vector<std::string>(), (niter - burnin) / thin );

cluster_indices_t initial_clusteringmvb (carcinoma.n_rows,1);

auto start_gibbs           = std::chrono::system_clock::now();
mixturemvb->fit(carcinoma , initial_clusteringmvb, false, priormvb , niter,  burnin,  thin, false, logger);
auto end_gibbs             = std::chrono::system_clock::now();
auto elapsed_gibbs         = end_gibbs - start_gibbs;
auto total_gibbs           = elapsed_gibbs.count() / 1000000.0;
COUT_STREAM << "Total time: " << total_gibbs << "ms" << std::endl ;

}


BOOST_AUTO_TEST_SUITE( test_suite_multivariate_binomial_dataset )

BOOST_AUTO_TEST_CASE( test_multivariate_binomial_100_10_10 )
{
	test_Mixture_MultivariateBinomial(100, 10, 10);
}


BOOST_AUTO_TEST_SUITE_END()




