/*
 * verbose.hpp
 *
 *  Created on: Jun 14, 2019
 *      Author: toky
 */

#ifndef ANTMAN_SRC_MIXTURE_CPP_VERBOSE_HPP_
#define ANTMAN_SRC_MIXTURE_CPP_VERBOSE_HPP_

#include <cstdlib>

#define EXTRA_LEVEL   4
#define DEBUG_LEVEL   3
#define INFO_LEVEL    2
#define LOG_LEVEL     1
#define WARNING_LEVEL 1
#define ERROR_LEVEL   0

extern int VERBOSE_LEVEL;

#ifdef Rcpp
#define COUT_STREAM Rcpp::Rcout
#define CERR_STREAM Rcpp::Rcerr

static inline void stop_cmd () {Rcpp::stop("Error inside the package.\n");}
#else
#define COUT_STREAM std::cout
#define CERR_STREAM std::cerr
static inline void stop_cmd () {abort();}
#endif

#ifdef VERBOSE_BINARY
#define VERBOSE_EXTRA(msg)                         {if (VERBOSE_LEVEL >= EXTRA_LEVEL)   COUT_STREAM  << msg << std::endl;};
#define VERBOSE_DEBUG(msg)                         {if (VERBOSE_LEVEL >= DEBUG_LEVEL)   COUT_STREAM  << msg << std::endl;};
#else
#define VERBOSE_EXTRA(msg)                         {};
#define VERBOSE_DEBUG(msg)                         {};
#endif

#define VERBOSE_INFO(msg)                          {if (VERBOSE_LEVEL >= INFO_LEVEL)    CERR_STREAM  << msg << std::endl;};
#define VERBOSE_LOG(msg)                           {if (VERBOSE_LEVEL >= LOG_LEVEL)     CERR_STREAM  << msg << std::endl;};
#define VERBOSE_WARNING(test,msg) {if (not (test)) {if (VERBOSE_LEVEL >= WARNING_LEVEL) CERR_STREAM  << msg << std::endl;}};
#define VERBOSE_ERROR(msg)                         {if (VERBOSE_LEVEL >= ERROR_LEVEL)   CERR_STREAM  << msg << std::endl;  stop_cmd () ;  };
#define VERBOSE_ASSERT(test,msg)  {if (not (test)) {                                    CERR_STREAM  << msg << std::endl;  stop_cmd () ; }};





#endif /* ANTMAN_SRC_MIXTURE_CPP_VERBOSE_HPP_ */
