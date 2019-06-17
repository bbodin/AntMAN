/*
 * maths.hpp
 *
 *  Created on: Jun 14, 2019
 *      Author: toky
 */

#ifndef ANTMAN_MIXTURE_CPP_MATHS_HPP_
#define ANTMAN_MIXTURE_CPP_MATHS_HPP_

#include "verbose.hpp"


#ifdef NO_RCPP

#include <armadillo>

#ifndef M_LN_2PI
#define M_LN_2PI	1.837877066409345483560659472811	/* log(2*pi) */
#endif

#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi))
								 == log(2*pi)/2 */
#endif

static inline double am_rpois   (double n)                        {
	 static std::default_random_engine generator;
	        std::poisson_distribution<int> distribution(n);
	return (distribution(generator));
}
static inline double am_runif   (double a ,double b)                   {
	 static std::default_random_engine generator;
	        std::uniform_real_distribution<double> distribution(a,b);
	 return (distribution(generator));
}

static inline double am_rnbinom (double a ,double b)                   {
	 static std::default_random_engine generator;
	 std::negative_binomial_distribution<int> distribution(a,b);
	 return (distribution(generator));
}

static inline double am_rgamma (double a ,double b)                   {
	 static std::default_random_engine generator;
	 std::gamma_distribution<double> distribution(a,b);
	 return (distribution(generator));
}
static inline double am_rnorm (double a ,double b)                   {
	 static std::default_random_engine generator;
	 std::normal_distribution<double> distribution(a,b);
	 return (distribution(generator));
}

static inline double am_rchisq (double a)                   {
	 static std::default_random_engine generator;
	 std::chi_squared_distribution<double> distribution(a);
	 return (distribution(generator));
}


static inline double am_pnorm   (double,double, double,bool,bool) {VERBOSE_ERROR("Unsupported function: am_pnorm  "); return 0.0;}
static inline double am_qnorm   (double,double, double,bool,bool) {VERBOSE_ERROR("Unsupported function: am_qnorm  "); return 0.0;}


#else
#include <RcppArmadillo.h>


static inline double am_rpois   (double n)                    {return R::rpois(n);  }
static inline double am_runif   (double a ,double b)          {return R::runif(a,b);  }
static inline double am_rnbinom (double a ,double b)          {return R::rnbinom(a,b);}
static inline double am_rgamma (double a ,double b)           {return R::rgamma(a,b); }
static inline double am_rnorm (double a ,double b)            {return R::rnorm(a,b); }
static inline double am_rchisq (double a)                     {return R::rchisq(a); }


static inline double am_pnorm   (double,double, double,bool,bool) {VERBOSE_ERROR("Unsupported function: am_pnorm  "); return 0.0;}
static inline double am_qnorm   (double,double, double,bool,bool) {VERBOSE_ERROR("Unsupported function: am_qnorm  "); return 0.0;}



#endif

#endif /* ANTMAN_MIXTURE_CPP_MATHS_HPP_ */