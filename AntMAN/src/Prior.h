/*
 * Prior.hpp
 *
 *  Created on: Apr 8, 2019
 */

#ifndef ANTMAN_SRC_PRIOR_H_
#define ANTMAN_SRC_PRIOR_H_

#include "math_utils.h"
#include "utils.h"
// --------------------------------------------------------------------------------------------------------------------


class h_param_t {
public:
#ifdef HAS_RCPP
	virtual const Rcpp::List get_Rcpp_list () const  = 0;
#endif
	virtual ~h_param_t() {};
};

class q_param_t {
public:
#ifdef HAS_RCPP
	virtual const Rcpp::List get_Rcpp_list () const  = 0;
#endif
	virtual ~q_param_t() {};
};

template<typename Q_t>
class gamma_h_param_t : public h_param_t {
public:
	bool gamma_is_fixed;
	double gamma;
	const double a,b; // hyper-prior parameters for h
	double lsd, lsd_g;       // this is the standard deviation of the MH algorithm to update gamma.
	gamma_h_param_t (double gamma, double a, double b, double lsd) : gamma_is_fixed (false), gamma (gamma) , a(a), b(b), lsd(lsd), lsd_g(1) {}
	gamma_h_param_t (              double a, double b, double lsd) : gamma_is_fixed (false), gamma (am_rgamma(a,b)) , a(a), b(b), lsd(lsd), lsd_g(1) {}
	gamma_h_param_t (double gamma) : gamma_is_fixed (true), gamma (gamma), a(0), b(0), lsd(0), lsd_g(1)  {}
#ifdef HAS_RCPP
	virtual const Rcpp::List get_Rcpp_list () const  {
		return Rcpp::List::create(Rcpp::Named("gamma") = this->gamma ) ;
	};
#endif


	void update (const  double U, const  int K, const std::vector<int> &nj , const  Q_t& q_param) {
		if (this->gamma_is_fixed) return;

		const double vecchio = this->gamma;
		const double lmedia = std::log(vecchio);

			//Propose a new value
			const double lnuovo=am_rnorm(lmedia,lsd);
			const double nuovo=std::exp(lnuovo);


	//		const double log_full_gamma_new = log_full_EPPF (nuovo , K , nj,  U , this->q_param.lambda ) + (ah-1)*std::log(nuovo)-bh*nuovo;
	//		const double log_full_gamma_vec = log_full_EPPF (vecchio , K , nj,  U , this->q_param.lambda ) + (ah-1)*std::log(vecchio)-bh*vecchio;
	//		const double ln_acp = (log_full_gamma_new - lmedia) - (log_full_gamma_vec - lnuovo);


			double ln_acp = (q_param.log_full_gamma(nuovo   , K , nj,   U , a, b) - lmedia)
					      - (q_param.log_full_gamma(vecchio , K , nj,   U , a, b) - lnuovo);


			const double lnu=std::log(am_runif(0.0,1.0));

			this->gamma = lnu<ln_acp ? nuovo : vecchio;

			lsd = update_lsd (  lsd,  ln_acp,  lsd_g++) ;

	}

};



class Prior {


public :
	virtual void    update (const double U, const int K, const std::vector<int> &nj ) = 0 ;

	virtual double  get_gamma() const = 0;
	virtual const h_param_t * get_h() const = 0;
	virtual const q_param_t * get_q() const = 0;

	virtual int     init_M_na(const int K)= 0;
	virtual int     update_M_na(const double U ,  const int K)= 0;

	virtual        ~Prior() {};

};

template <typename H_t, typename Q_t>
class TypedPrior : public Prior {

protected :
	H_t h_param;
	Q_t q_param;
public:


	void update (const double U, const int K, const std::vector<int> &nj ) {
		 this->q_param.update (U, K, this->h_param);
		 this->h_param.update (U, K, nj, this->q_param);
	};

	double get_gamma() const {return this->h_param.gamma;};

	const h_param_t* get_h() const {return & (this->h_param);};
	const q_param_t* get_q() const {return & (this->q_param);};

	TypedPrior(H_t h_param, Q_t q_param)             : h_param(h_param), q_param(q_param) {};


	virtual    ~TypedPrior() {};

};


#endif /* ANTMAN_SRC_PRIOR_H_ */
