/*
 * MixtureUnivariateBernoulli.hpp
 *
 *  Created on: Mar 8, 2019
 *      Author: toky
 */

#ifndef PROBITFMMNEW_SRC_MIXTUREUNIVARIATEBERNOULLI_HPP_
#define PROBITFMMNEW_SRC_MIXTUREUNIVARIATEBERNOULLI_HPP_


#include <RcppArmadillo.h>
#include "Mixture.hpp"

class Mixture_UnivariateBernoulli : public UnivariateMixture  {

	//ParametricPrior
	int         _mb;
	double _a0, _b0;

	//Tau
	std::vector <double> _theta;

public:
	typedef arma::vec input_t;

public :
	Mixture_UnivariateBernoulli (const double a0, const double b0, const int mb) :  _mb(mb), _a0 (a0), _b0 (b0) {}

	virtual void init_tau (const input_t & y, const int M) {
		 _theta .resize(M)  ;

		 const double b0 = _b0;
		 const double a0 = _a0;
		 const double mb = _mb;

				for(int l=0;l<M;l++){
					_theta[l] = R::rbeta(a0, b0);
				}

	}

		virtual cluster_indices_t  up_ci(const  input_t & y,
				const long M,
				const Rcpp::NumericVector & S_current) {

			const int    mb    = _mb;
			const int    n     = y.size(); // TODO[DISCUSS ME]: there is a problem between colvec rowvec and matrix and the way data is layout (d,n) or (n,1)

			const std::vector<double>& theta_current = _theta;
			Rcpp::NumericVector Log_S_current = log(S_current);
			cluster_indices_t ci_current(n);
			Rcpp::NumericVector random_u   = Rcpp::runif(n,0.0,1.0 );

			for (int i=0; i < n; i++) {

				arma::vec pesi(M);
				double max_lpesi=-INFINITY;

				for(int l=0;l<M;l++){

					 double ldensi =   y[i] * log (theta_current[l])  + (mb - y[i]) * log (1 - theta_current[l]);

					 pesi[l]=Log_S_current[l]+ldensi;
					 if(max_lpesi<pesi[l]){max_lpesi=pesi[l];}
				}

				// I put the weights in natural scale and re-normalize then
				pesi = exp(pesi - max_lpesi);
				pesi = pesi / sum(pesi);

				const double u = random_u[i];
				double cdf = 0.0;
				unsigned int ii = 0;
				while (u >= cdf) { // This loop assumes (correctly) that R::runif(0,1) never return 1.
					cdf += pesi[ii++];
				}
				ci_current[i] = ii;
			}

			return  ci_current ;



		}

		 virtual allocation_result up_allocated_nonallocated (
				const int K ,
				const int M ,
				const cluster_indices_t & ci_current ,
				const cluster_indices_t & ci_star  ,
				const double gamma_current,
				const double U_current,
				const  input_t & y ) {

				const long   n     = y.size();
				const double a0    = _a0;
				const double b0    = _b0;
				const int    mb    = _mb;


				//Allocation_result output ;
				std::vector<double> theta_current(M);
				Rcpp::NumericVector S_current    = Rcpp::NumericVector(M);

				cluster_indices_t ci_reorder(y.size());
				ci_reorder.fill(-1);
				std::vector<int>    nj(K);
				std::map< int, std::vector<int> > clusters_indices;

				for(int i=0;i<n;i++){
					clusters_indices[ci_current[i]].push_back(i);
				}

				for(int local_index=0;local_index<K;local_index++){
					const int key = ci_star[local_index];
					nj[local_index] = clusters_indices[key].size();
					for (auto v : clusters_indices[key]) {
						ci_reorder[v]=local_index;
					}
				}

				for(int l=0; l<K;l++){
					// Find the index of the data in the l-th cluster
					std::vector<int> & which_ind=clusters_indices [ci_star[l]];

					//Prepare the variable that will contain the data in the cluster
					arma::colvec y_l = arma::colvec(nj[l]);
					//Separate the data in each cluster and rename the cluster
					for(int it=0;it<nj[l];it++){
						y_l[it]=y[which_ind[it]];

					}

					// Call the parametric function that update the
					// parameter in each cluster
					// See the boxes 4.c.i and 4.b.i of
					// Figure 1 in the paper

					const int njl = y_l.n_elem; // This is the number of data in the cluster. I hope so


					//Since in our case the full conditionals are in closed form
					//they are Normal-inverse-gamma



					//First compute the posterior parameters

					const double ysum= arma::sum(y_l);

					//double ysum= arma::accu(y); //is there any difference between sun and accu un Armadillo?


					const double an = a0 + ysum;
					const double bn = njl  * mb - ysum + b0;
					theta_current[l] = R::rbeta (an, bn);

					// TODO[OPTIMIZE ME] : Cannot split or the random value are completly different
					// Update the Jumps of the allocated part of the process
					S_current[l]=R::rgamma(nj[l]+gamma_current,1./(U_current+1.0));

				}

				// Fill non-allocated

				for(int l=K; l<M;l++){
					theta_current[l] = R::rbeta (a0, b0) ;
					// TODO[CHECK ME] : theta_current must be non-zero
					S_current[l]=R::rgamma(gamma_current,1./(U_current+1.0));
				}



				_theta = theta_current;
				return allocation_result(ci_reorder , nj , S_current);




		 }
};




#endif /* PROBITFMMNEW_SRC_MIXTUREUNIVARIATEBERNOULLI_HPP_ */
