/*
 * utils.hpp
 *
 *  Created on: Jan 30, 2019
 *      Author: toky
 */

#ifndef PROBITFMM_SRC_UTILS_HPP_
#define PROBITFMM_SRC_UTILS_HPP_

#define DEBUG_LEVEL   4
#define EXTRA_LEVEL   3
#define INFO_LEVEL    2
#define WARNING_LEVEL 1

extern int VERBOSE_LEVEL;

#define VERBOSE_DEBUG(msg) {if (VERBOSE_LEVEL > DEBUG_LEVEL)  Rcpp::Rcout  << msg << std::endl;};
#define VERBOSE_INFO(msg)  {if (VERBOSE_LEVEL > INFO_LEVEL) Rcpp::Rcerr  << msg << std::endl;};
#define VERBOSE_ERROR(msg) {Rcpp::Rcerr << msg << std::endl;  Rcpp::stop("Error inside the package.\n"); };
#define VERBOSE_WARNING(test,msg) {if (not (test)) { if (VERBOSE_LEVEL >= WARNING_LEVEL) Rcpp::Rcerr << msg << std::endl;  }};
#define VERBOSE_ASSERT(test,msg) {if (not (test)) { Rcpp::Rcerr << msg << std::endl;  Rcpp::stop("Error inside the package.\n"); }};




typedef arma::ivec cluster_indices_t;

/// this is a function to convert an arma vec to a NumericVector object
template <typename T>
inline Rcpp::NumericVector arma2vec(const T& x) {
	return Rcpp::NumericVector(x.begin(), x.end());
}

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
	int ncols = Sig.n_cols;
	arma::vec Y = Rcpp::as<arma::vec>(Rcpp::rnorm(ncols,0.0,1.0 ));
	return mu +  arma::chol(Sig) * Y;
}

// TODO[LICENCE ISSUE !!] : We took that from Rcpp-dist, our code must be GPL !

inline arma::vec dmvnorm(const arma::mat& x, const arma::vec& mu,
        const arma::mat& S, const bool log_p = false) {

    arma::uword n = x.n_rows, m = x.n_cols;
    double det_S = arma::det(S);

    arma::mat S_inv = S.i();

    arma::vec result(n);

    arma::rowvec X(m);

    arma::rowvec Mu = mu.t();

    if ( log_p ) {
        double P = -1.0 * (x.n_cols/2.0) * M_LN_2PI - 0.5 * log(det_S);
        for ( arma::uword i = 0; i < n; ++i ) {

            X = x.row(i) - Mu;

            result[i] = arma::as_scalar(P - 0.5 * X * S_inv * X.t());
        }
        return result;
    }
    double P = 1.0 / sqrt(pow(M_2PI, m) * det_S);
    for ( arma::uword i = 0; i < n; ++i ) {
        X = x.row(i) - Mu;
        result[i] = arma::as_scalar(P * exp(-0.5 * X * S_inv * X.t()));
    }
    return result;
}

inline double dmvnorm1(const arma::vec& x, const arma::vec& mu,
        const arma::mat& S, const bool log_p = false) {

	return dmvnorm(x.t(), mu,S,log_p)[0];
}

inline double dmvnorm_raffa(const arma::vec& x, const arma::vec& mu,
        const arma::mat& S, const bool log_p = false) {
    return 0.0;

}


// TODO[LICENCE ISSUE !!] : We took that from Rcpp-dist, our code must be GPL !
inline arma::mat rwish(const int df, const arma::mat& S) {
    arma::uword m = S.n_cols;
    arma::uword i, j;
    arma::mat A(m, m, arma::fill::zeros);
    for ( i = 1; i < m; ++i ) {
        for ( j = 0; j < i; ++j ) {
            A(i, j) = R::rnorm(0.0, 1.0);
        }
    }
    for ( i = 0; i < m; ++i ) {
        A(i, i) = sqrt(R::rchisq(df - i));
    }
    arma::mat B = A.t() * arma::chol(S);
    return B.t() * B;
}

// TODO[LICENCE ISSUE !!] : We took that from Rcpp-dist, our code must be GPL !
inline arma::mat riwish(const int df, const arma::mat& S) {
    return rwish(df, S.i()).i();
}


//Proposition 2.3 og Robert
inline double rnorm_soprasoglia_cr2(double mumeno, double alpha=0){

	if(alpha==0){
		alpha=mumeno;
	}

	if(mumeno==0){
		Rcpp::Rcout<<__LINE__<<"Errore! mumeno non puo' essere zero\n";
	}

	while(true){

		// Step 1
		const double u1= R::runif(0,1);
		const double z = mumeno-std::log(u1)/alpha;

		//Step 2
		const double lM = std::pow(alpha,2)/2;
		const double laccept = -std::pow(z,2)/2+alpha*z-lM;

		// Step 3
		// TODO[CHECK ME] : could never finish ... ?
		const double u = R::runif(0,1);
		const double log_u = std::log(u);

		if(log_u<laccept){
			return z;
		}
	}
}



inline double fast_rnorm_truncated(const double mean, const bool sopra){


	double out;
	int factor = sopra ? 1 : -1;
	double imean = mean*factor;
	double soglia = -imean;

	if(soglia<0.1){
		double down = R::pnorm(soglia,0.0,1.0,true,false);
		double u = R::runif(down,1);
		out =R::qnorm(u, 0.0, 1.0,true,false);
	}
	else{
		out=rnorm_soprasoglia_cr2(soglia,0);
	}

	return( factor*(imean+out));
}


inline double rnorm_truncated(const double mean,const double sd, const double soglia, const bool sopra) {

	double out;

	int factor = sopra ? 1 : -1;
	const double meanf = mean*factor;

	const double tmp_soglia = (soglia*factor-meanf)/sd;

	if(tmp_soglia<0.1){
		double down = R::pnorm(tmp_soglia,0.0,1.0,true,false);
		double u = R::runif(down,1);
		out =R::qnorm(u, 0.0, 1.0,true,false);
	}
	else{
		out=rnorm_soprasoglia_cr2(tmp_soglia,0);
	}

	return( factor*(meanf+sd*out));
}





inline std::vector<int> which_eq(std::vector<int> x, int j) {
	int nx = x.size();
	std::vector<int> out;
	out.reserve(nx);
	for(int i = 0; i < nx; i++) {
		if (x[i]==j ) out.push_back(i);
	}
	return out;
}



////// Here I consider my own functions to sample from a discrete!

inline arma::vec sample_raf(unsigned int max,int hm, arma::vec weights, int plus1 ) {
	//this function sample hm observations from a discrete
	//distribution with support in 0,1,...,max-1
	//with p.m.f. proportional to weights
	//  if plus1 = 1 -> support is (1, max)

	double somma;
	double cdf;
	arma::vec out(hm);
	somma = sum(weights);


	if(!(somma>0)){
		Rcpp::Rcout<<"Problem with the sum of the weights!"<<"\n";
		throw std::runtime_error(" ");
	}

	if(max!=weights.size())
	{
		Rcpp::Rcout<<"The support has a length different from the weights"<<"\n";
		throw std::runtime_error(" ");
	}

	for(int g = 0; g < hm; g++){

		const double u = R::runif(0,1);

		cdf = 0.0;
		for(unsigned int ii = 0; ii < max; ii++){
			cdf += weights[ii]/somma;
			if(u < cdf){
				out[g] = ii + plus1;
				break;
			}
		}
	}

	return out;
}

inline double fast_sample_raf(const arma::vec& weights ) { // sample_raf when (plus1 == 0) and (hm == 1)
	const double somma = sum(weights); // TODO[CHECK ME] : Must be 1 not ??

	if(!(somma>0)){
		Rcpp::Rcout<<"Problem with the sum of the weights!"<<"\n";
		throw std::runtime_error(" ");
	}



	const double u = R::runif(0,1);

	double cdf = 0.0;
	for(unsigned int ii = 0; ii < weights.size(); ii++){
		cdf += weights[ii]/somma;
		if(u < cdf){
			return ii;
		}
	}
	// TODO[CHECK ME] What happen if condition not met ?
	VERBOSE_ERROR("This might be a mistake in the algorithm implementation.");
	return 0.0;
}





///////////////////////////////////////////////////////////////////////
//THE FOLLOWING ARE TWO FUNCTIONS TO SAMPLE FROM TRUNCATED GAUSSIAN////
///////////////////////////////////////////////////////////////////////






#endif /* PROBITFMM_SRC_UTILS_HPP_ */
