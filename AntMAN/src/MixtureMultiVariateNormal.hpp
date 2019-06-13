/*
 * MixtureMultiVariateNormal.hpp
 *
 *  Created on: Mar 8, 2019
 *      Author: toky
 */

#ifndef PROBITFMMNEW_SRC_MIXTUREMULTIVARIATENORMAL_HPP_
#define PROBITFMMNEW_SRC_MIXTUREMULTIVARIATENORMAL_HPP_


#include <RcppArmadillo.h>
#include "utils.hpp"
#include "Mixture.hpp"
// [[Rcpp::depends(RcppArmadillo)]]

class Mixture_MultivariateNormal: public MultivariateMixture  {

	// ParametricPrior
	arma::vec    _mu0;
	double       _ka0;
	unsigned int _nu0;
	arma::mat    _Lam0;

	// Tau
	arma::mat  _mu_current;
	arma::cube _Sig_current;


public :
	Mixture_MultivariateNormal (const arma::vec & mu0, const double ka0, const unsigned int nu0, const arma::mat & Lam0) : _mu0 (mu0), _ka0 (ka0), _nu0 (nu0), _Lam0  (Lam0){}

	virtual void init_tau (const input_t & y, const int M){

		VERBOSE_DEBUG(" init_tau (const input_t & y, const int M)");

		arma::mat   mu_current  (M,y.n_cols);
		arma::cube  Sig_current (y.n_cols,y.n_cols, M) ;

		 _mu_current = mu_current ;
		 _Sig_current= Sig_current ;

				const arma::colvec mu0  = _mu0;
				const double       ka0  = _ka0;
				const unsigned int nu0  = _nu0;
				const arma::mat    Lam0 = _Lam0;


				for(int l=0;l<M;l++){

					//  TODO[CHECK ME] : arma::mat rwish(const int df, const arma::mat& S)
					//  arma::vec mvrnormArma(arma::colvec mu, arma::mat Sig)
					// TODO[CHECK ME] : Raffa study the lam0 - 1
					const arma::mat res = riwish (nu0, Lam0);

					_Sig_current.slice(l) = res;

					const arma::colvec tmp = mvrnormArma (mu0, res / ka0) ;

					_mu_current.row(l)   =  tmp.t();

				}

				VERBOSE_DEBUG(" init_tau finished");

		 }

		virtual cluster_indices_t  up_ci(const  input_t & y,
				const long M,
				const Rcpp::NumericVector & S_current){


			VERBOSE_DEBUG("GibbsFramework<Tau_MultivariateNormal>::up_ci");

			const int n = y.n_rows;

			const arma::mat& mu_current = _mu_current;
			const arma::cube& Sig_current = _Sig_current;

			cluster_indices_t ci_current(n);
			Rcpp::NumericVector Log_S_current = log(S_current);
			Rcpp::NumericVector random_u   = Rcpp::runif(n,0.0,1.0 );

			VERBOSE_DEBUG("GibbsFramework<Tau_MultivariateNormal>::up_ci: for (int i=0; i < n; i++) {");

			VERBOSE_DEBUG("y  :" << y.n_rows << "x" << y.n_cols);
			VERBOSE_DEBUG("mu_current :" << mu_current.n_rows << "x" << mu_current.n_cols);
			VERBOSE_DEBUG("Sig_current  :"  << Sig_current.n_rows << "x" << Sig_current.n_cols << "x" << Sig_current.n_slices );


			omp_set_dynamic(0);

			#pragma omp parallel for num_threads(8)
			for (int i=0; i < n; i++) {

				arma::vec pesi(M);
				double max_lpesi=-INFINITY;

				for(int l=0;l<M;l++){

					const arma::vec dmvnorm1_mu  = mu_current.row(l).t();
					const arma::mat dmvnorm1_Sig =Sig_current.slice(l);


					double ldensi = dmvnorm(y.row(i), dmvnorm1_mu,dmvnorm1_Sig,true)[0];
					 pesi[l]=Log_S_current[l]+ldensi;
					 if(max_lpesi<pesi[l]) {
						 max_lpesi=pesi[l];
					 }
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


				const int n = y.n_rows;
				const int d = y.n_cols;

				const arma::vec    mu0  = _mu0;
				const double       ka0  = _ka0;
				const unsigned int nu0  = _nu0;
				const arma::mat    Lam0 = _Lam0;


				 arma::mat mu_current (M, d);
				 arma::cube Sig_current (d , d, M);

				Rcpp::NumericVector S_current    = Rcpp::NumericVector(M);

				cluster_indices_t ci_reorder(y.n_rows);
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

				VERBOSE_DEBUG("I know but does it print ? ");
				omp_set_dynamic(0);

				for(int l=0; l<K;l++){

					VERBOSE_DEBUG("++");
					// Find the index of the data in the l-th cluster
					std::vector<int> & which_ind=clusters_indices [ci_star[l]];

					VERBOSE_DEBUG("++");
					//Prepare the variable that will contain the data in the cluster
					std::vector <arma::vec>  y_l (nj[l]);
					//Separate the data in each cluster and rename the cluster
					for(int it=0;it<nj[l];it++){
						y_l[it]=y.row(which_ind[it]).t();

					}

					// Call the parametric function that update the
					// parameter in each cluster
					// See the boxes 4.c.i and 4.b.i of
					// Figure 1 in the paper

					VERBOSE_DEBUG("++");
					const int njl = y_l.size(); // This is the number of data in the cluster. I hope so

					VERBOSE_DEBUG("++");

					//Since in our case the full conditionals are in closed form
					//they are Normal-inverse-gamma


					//Firs compute the posterior parameters

					VERBOSE_DEBUG("y_l is a std::vector of  =>" << y_l.size());

					const arma::vec ysum = y_l.size() ? vectorsum(y_l) : arma::zeros(d);

					VERBOSE_DEBUG("ysum =>" << ysum.n_rows << "x" << ysum.n_cols);
					VERBOSE_DEBUG("++");
					const double ka0nokon = (ka0 * njl ) / (ka0 + njl);

					const arma::vec ybar =(ysum / (double) njl); // cast?

					VERBOSE_DEBUG("ybar =>" << ybar.n_rows << "x" << ybar.n_cols);
					VERBOSE_DEBUG("mu0  =>" << mu0.n_rows << "x" << mu0.n_cols);
					const arma::vec ybarminusmu0 = (ybar - mu0);
					const arma::rowvec ybarminusmu0t = ybarminusmu0.t();
					const arma::mat  ybarmatmul = ybarminusmu0 * ybarminusmu0t;

					VERBOSE_DEBUG("++");
					arma::mat S2  = arma::zeros(d, d);
					VERBOSE_DEBUG("S2 = " << S2);
					for (int i = 0 ; i < njl ; i ++) {
						const arma::vec ylmyb = (y_l[i] - ybar);


						const arma::rowvec ylmybt = ylmyb.t();

						const arma::mat matmul = ylmyb * ylmybt;
						VERBOSE_DEBUG("matmul = " << matmul);
						VERBOSE_DEBUG("S2 = " << S2);
						S2 =  S2 + matmul; // TODO[CHECK ME] : he does not trust me :-(
					}

					// Then the parameters of the posteriorr  Normal-inverse-gamma a posteriori are
					const double kan = ka0  +  (double) njl; // maybe the cast at double in the sum is not needed
					const unsigned int nun = nu0 + njl;  // maybe the cast at  double in the sum is not needed
					const arma::vec mun = (ysum+mu0*ka0)/kan;
					const arma::mat Lamn = Lam0 + S2 + ka0nokon * ybarmatmul;

					const arma::mat Sig_l =  riwish (nun, Lamn) ;
					const arma::vec mu_l   = mvrnormArma(mun, Sig_l/kan);


					Sig_current.slice(l) = Sig_l; // In case 4 we have to update a matrix
					mu_current.row(l)  = mu_l.t();
				}

				for(int l=0; l<K;l++){
					// TODO[CHECK ME] : I split the loop, random generation is not the saame as before, but it should be theoricaly equivalent, isnt it ?
					// Update the Jumps of the allocated part of the process
					S_current[l] = R::rgamma(nj[l]+gamma_current,1./(U_current+1.0));

				}

				VERBOSE_DEBUG("--");
				// Fill non-allocated

				for(int l=K; l<M;l++){

							const arma::mat res = riwish (nu0, Lam0);
							Sig_current.slice(l) = res; // TODO[CHECK ME] : Raffa study the lam0 - 1
							mu_current.row(l)   = mvrnormArma (mu0, Sig_current.slice(l) / ka0).t() ;
							S_current[l]=R::rgamma(gamma_current,1./(U_current+1.0));
						}


			 // TODO[OPTIMIZE ME] : same as initializatrion but initialization can be anything philosophical problem !!!


				this->_Sig_current = Sig_current;
				this->_mu_current = mu_current;

				return allocation_result(ci_reorder , nj , S_current);



		 }
};






#endif /* PROBITFMMNEW_SRC_MIXTUREMULTIVARIATENORMAL_HPP_ */
