/*
 * MixtureUnivariateProbit.hpp
 *
 *  Created on: Mar 8, 2019
 *      Author: toky
 */

#ifndef PROBITFMMNEW_SRC_MIXTUREUNIVARIATEPROBIT_HPP_
#define PROBITFMMNEW_SRC_MIXTUREUNIVARIATEPROBIT_HPP_


#include <RcppArmadillo.h>
#include "utils.hpp"
#include "Mixture.hpp"
#include <omp.h>

// This is the parametric update on the data of cluster l
static inline arma::colvec up_beta_parametric(
		const int hmany ,
		const int k_x,
		const arma::ucolvec & y_l ,
		const  arma::mat    & Sig_betam1,
		const  arma::mat    & X_l,
		const arma::colvec  & mu_beta,
		const arma::colvec  & old_z_l){

	arma::colvec z_l = old_z_l;

	arma::colvec beta_loc;
	int njl = y_l.n_elem;

	const arma::mat tX_l     = arma::trans(X_l);
	const arma::mat Vbeta_l  = arma::pinv(Sig_betam1 + tX_l*X_l);
	const arma::mat sigmu    = Sig_betam1 * mu_beta ;
	const int ncols          = Vbeta_l.n_cols;
	const arma::mat CVbeta_l = arma::chol(Vbeta_l);

	for(int hhm=0; hhm < hmany; hhm++){

		/// Mean of the full conditional of the beta
		const arma::mat Ebeta_l   = Vbeta_l * (sigmu + tX_l * z_l);
		const arma::vec Y         = Rcpp::as<arma::vec>(Rcpp::rnorm(ncols,0.0,1.0 ));
		beta_loc =  Ebeta_l +  CVbeta_l * Y;
		const arma::colvec matmut = X_l * beta_loc;

		for(int ijl=0;ijl<njl;ijl++){
			z_l[ijl]= rnorm_truncated(matmut[ijl],1,0,y_l[ijl]==1);
		}// Close the loop on the z_loc

	}/// Close the loop on how many

	return beta_loc;

}

class Mixture_UnivariateProbit: public UnivariateMixture {

	// Parametric Prior
	arma::mat     _X;
	arma::colvec  _mu_beta;
	arma::mat     _Sig_beta;

	//Tau
	arma::mat           _Beta_current;
	arma::colvec        _z_current ;

	// cached values
	arma::mat     _Sig_betam1;
	int           _k_x;

public :
	Mixture_UnivariateProbit (const arma::mat & X, const arma::colvec & Mu, const arma::mat & Sig) : _X(X), _mu_beta(Mu), _Sig_beta(Sig) , _Sig_betam1 ( arma::pinv(Sig)),  _k_x  (X.n_cols){}

	virtual void init_tau (const input_t & y, const int M) {
		_Beta_current.resize (M,_k_x);
		_z_current.resize (y.size()) ;

		for(int it=0;it<M;it++){
			_Beta_current.row(it)= mvrnormArma(_mu_beta,_Sig_beta);
		}

		for(unsigned int i=0;i<y.size();i++){
			_z_current[i]=rnorm_truncated(0,1,0,(y[i]==1)?1:0);
		}

	}

		virtual cluster_indices_t  up_ci(const  input_t & y,
				const long M,
				const Rcpp::NumericVector & S_current) {


			const int n = y.size();

			const arma::mat    Beta_current = _Beta_current;
			const arma::colvec z_current = _z_current;
			const arma::mat  X = _X;


			cluster_indices_t ci_current(n);

			Rcpp::NumericVector Log_S_current = log(S_current);

			static double const SQRT_2PI = 0.918938533204672741780329736406 ;
			arma::mat matrixProduct = X*Beta_current.t();

			std::vector<arma::vec> all_pesi(n);

			// A classic "#pragma omp parallel for" enable speedup but is unstable.
			// With this code, OpenMP's speed up is stable and satisfactory (on my laptop...).
			omp_set_dynamic(0);
		    #pragma omp parallel for num_threads(4)
			for(int i=0;i<n;i++){
				arma::vec pesi(M);
				double max_lpesi=-INFINITY;

				for(int l=0;l<M;l++){
					double linpred_local= matrixProduct(i,l) ;

					// TODO : Parametric here ! ldensi change.
					double x     = fabs ((z_current[i] - linpred_local));
					double ldensi=	-(SQRT_2PI + 0.5 * x * x );
					// Speed-up 30% using homemade dnorm
					//double ldensi2=R::dnorm(z_current[i],linpred_local,1,true);

					pesi[l]=Log_S_current[l]+ldensi;
					if(max_lpesi<pesi[l]){max_lpesi=pesi[l];}
				}

				// I put the weights in natural scale and re-normalize then
				pesi = exp(pesi - max_lpesi);
				pesi = pesi / sum(pesi);
				all_pesi[i] = pesi;
			}

			// TODO: Does this part needs to be in order because of the rnorm ?
			for(int i=0;i<n;i++){
				//Let me sample the allocation of the i-th individual
				ci_current[i]= fast_sample_raf(all_pesi[i]);
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


				const long n = y.size();
				const long k_x = _k_x;
				const arma::mat  X = _X;
				const arma::colvec z_current = _z_current;
				const arma::mat  Sig_betam1 = _Sig_betam1;
				const arma::colvec mu_beta = _mu_beta;


				const arma::mat  Sig_beta = _Sig_beta;


				//Allocation_result output ;
				arma::mat           Beta_current (M,k_x);
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
					arma::mat X_l = arma::mat(nj[l],k_x);
					arma::colvec z_l = arma::colvec(nj[l]);
					arma::ucolvec y_l = arma::ucolvec(nj[l]);

					//TODO : This shuffling is not required
					//Separate the data in each cluster and rename the cluster
					for(int it=0;it<nj[l];it++){
						X_l.row(it)=X.row(which_ind[it]);
						z_l[it]=z_current[which_ind[it]];
						y_l[it]=y[which_ind[it]];

					}

					// Call the parametric function that update the
					// parameter in each cluster
					// See the boxes 4.c.i and 4.b.i of
					// Figure 1 in the paper

					arma::colvec beta_l= up_beta_parametric(10 , k_x, y_l ,   Sig_betam1, X_l, mu_beta,  z_l);

					Beta_current.row(l)= beta_l;


					// TODO : Cannot split or the random value are completly different
					// Update the Jumps of the allocated part of the process
					S_current[l]=R::rgamma(nj[l]+gamma_current,1./(U_current+1.0));

				}

				// Fill non-allocated


				for(int l=K; l<M;l++){
					arma::colvec beta_l = mvrnormArma(mu_beta,Sig_beta);
					Beta_current.row(l)= beta_l;
					S_current[l]=R::rgamma(gamma_current,1./(U_current+1));
				}


				arma::colvec new_z_current(n);
				for(int i=0;i<n;i++){
					arma::colvec beta_l= Beta_current.row(ci_reorder[i]);
					double v = arma::as_scalar(X.row(i)*beta_l);
					new_z_current[i]=rnorm_truncated(v,1,0,y[i]==1);
				}


				_Beta_current = Beta_current;
				_z_current = new_z_current;

				return allocation_result (ci_reorder , nj , S_current);



		 }



};

#endif /* PROBITFMMNEW_SRC_MIXTUREUNIVARIATEPROBIT_HPP_ */
