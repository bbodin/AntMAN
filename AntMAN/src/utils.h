/*
 * utils.hpp
 *
 *  Created on: Jan 30, 2019
 */

#ifndef ANTMAN_SRC_UTILS_HPP_
#define ANTMAN_SRC_UTILS_HPP_

#include "math_utils.h"
#include "verbose.h"

typedef arma::ivec cluster_indices_t;

inline arma::vec vectorsum (std::vector <arma::vec> elems ) {
	arma::vec out = elems[0];
	for (unsigned int i = 1 ; i < elems.size() ; i ++) {
		out += elems[i] ;
	}
	return out;
}

inline double  update_lsd ( double lsd, double ln_acp, double iter) {


	// This is a new parameter to adjust lsd (ADAPTIVE METROPOLIS; the user could be allowed to set a different value in (-1,0) different that -0.7; Even if it is dangerous to change it
	double wg=std::pow(iter,-0.7);
	// This is a new parameter to adjust lsd (ADAPTIVE METROPOLIS; the user could be allowed to set a different value in (0,1) different that 0.234; We should worn however.
	// The adaptive Rejection Metropolis Hasting we are going to use is Algorithm 5 of Griffin Stephens (2003)
	double bartau =0.234;


	// Adaptive metropolis : Algorithm 5 Griffin Sthephens
	//lsd= lsd+wg*(std::exp(std::min(0.0,ln_acp))-bartau);
	/* Andrea suggestions Apparently Griffin uses the multiplicative version
	Per le'algoritmo 5 e parte del 6, anziche' scrivere:
	s2(new) = s2(old) + wg*(alpha - tau)
	scrivi:
	s2(new) = s2(old) * exp( wg*(alpha - tau) )
	*/
	lsd = lsd * std::exp(wg*(std::exp(std::min(0.0,ln_acp))-bartau));

	if(lsd<std::pow(10,-50)){
		lsd=std::pow(10,-50);
	}
	if(lsd>std::pow(10,50)){
		lsd=std::pow(10,50);
	}


	return lsd;

}


/// This is a function to sample one observation from a multivariate
//  Normal with mean vector mu and varcov Sig

inline arma::vec mvrnormArma(arma::colvec mu, arma::mat Sig) {

	VERBOSE_ASSERT(Sig.is_sympd(), "mvrnormArma requires Sig to be symmetric. It is not S = " << Sig);

	arma::vec Y = arma::randn<arma::vec>(Sig.n_cols);

	return mu +  arma::chol(Sig) * Y;
}


inline double dmvnormZero(const arma::mat& x, const arma::vec& mu, const arma::mat& S, const bool log_p = false) {

    arma::uword m = x.n_cols;
    double S_det = arma::det(S);
    arma::mat S_inv = arma::inv(S);
    arma::rowvec X(m);
    arma::rowvec Mu = mu.t();
    X = x.row(0) - Mu;
    if ( log_p ) {
        double P = -1.0 * (x.n_cols/2.0) * M_LN_2PI - 0.5 * log(S_det);
        return arma::as_scalar(P - 0.5 * X * S_inv * X.t());
    } else {
    	double P = 1.0 / sqrt(pow(M_2PI, x.n_cols) * S_det);
    	return arma::as_scalar(P * exp(-0.5 * X * S_inv * X.t()));
    }

}

// Strongly inspired from RcppDist, but assume this does not impose the GPL.
inline arma::mat riwish(const int df, const arma::mat& iS) {

	arma::mat S = arma::inv(iS);

	VERBOSE_ASSERT(S.is_sympd(), "riwish requires S to be symmetric. It is not S = " << S << " and iS = " << S);

	arma::uword m = S.n_cols;

    arma::mat A(m, m, arma::fill::zeros);

    for (arma::uword i = 1; i < m; ++i ) {
    	//A.col(i) = Rcpp::as<arma::vec>(Rcpp::rnorm(i)); // Need to test that
    	for (arma::uword j = 0; j < i; ++j ) {
    		A(i, j) =  am_rnorm(0.0, 1.0);
    	}
    }
    for (arma::uword i = 0; i < m; ++i ) {
    	A(i, i) = sqrt(am_rchisq(df - i));
    }
    arma::mat B = A.t() * arma::chol(S);

    return  arma::inv(B.t() * B);
}




#endif /* ANTMAN_SRC_UTILS_HPP_ */
