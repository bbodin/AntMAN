/*
 * verbose.hpp
 *
 *  Created on: Jun 14, 2019
 */

#ifndef ANTMAN_SRC_MIXTURE_CPP_VERBOSE_HPP_
#define ANTMAN_SRC_MIXTURE_CPP_VERBOSE_HPP_


#ifdef HAS_RCPP
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#define COUT_STREAM Rcpp::Rcout
#define CERR_STREAM Rcpp::Rcerr
static inline void stop_cmd () {Rcpp::stop("Error inside the package.\n");}
#else
#include <cstdlib>
#define COUT_STREAM std::cout
#define CERR_STREAM std::cerr
static inline void stop_cmd () {abort();} // Commented for R Package
#endif

#define VERBOSE_COLOR true

#define PURPLE_COLOR (VERBOSE_COLOR?"\033[1;35m":"")
#define RED_COLOR    (VERBOSE_COLOR?"\033[1;31m":"")
#define YELLOW_COLOR (VERBOSE_COLOR?"\033[0;33m":"")
#define GREEN_COLOR  (VERBOSE_COLOR?"\033[1;32m":"")
#define BLUE_COLOR   (VERBOSE_COLOR?"\033[1;34m":"")
#define RESET_COLOR  (VERBOSE_COLOR?"\033[0m":"")

#define EXTRA_LEVEL   4
#define DEBUG_LEVEL   3
#define INFO_LEVEL    2
#define LOG_LEVEL     1
#define WARNING_LEVEL 1
#define ERROR_LEVEL   0

extern int VERBOSE_LEVEL;

#define VERBOSE_GENERIC_MSG(thr, out, color, msg)      {if (VERBOSE_LEVEL >= thr)    out  << "[" << thr << "] " << color  << msg << RESET_COLOR << std::endl;              };
#define VERBOSE_GENERIC_END(thr, out, color, msg)      {if (VERBOSE_LEVEL >= thr)    out  << "[" << thr << "] " << color  << msg << RESET_COLOR << std::endl; stop_cmd () ;};


#ifdef VERBOSE_BINARY
#define VERBOSE_EXTRA(msg)                          VERBOSE_GENERIC_MSG(EXTRA_LEVEL,    CERR_STREAM, BLUE_COLOR,  msg)
#define VERBOSE_DEBUG(msg)                          VERBOSE_GENERIC_MSG(DEBUG_LEVEL,    CERR_STREAM, BLUE_COLOR,  msg)
#else
#define VERBOSE_EXTRA(msg)                         {};
#define VERBOSE_DEBUG(msg)                         {};
#endif

#define VERBOSE_INFO(msg)                            VERBOSE_GENERIC_MSG(INFO_LEVEL,    CERR_STREAM, GREEN_COLOR,  msg)
#define VERBOSE_LOG(msg)                             VERBOSE_GENERIC_MSG(LOG_LEVEL,     CERR_STREAM, RESET_COLOR,  msg)
#define VERBOSE_WARNING(msg)                         VERBOSE_GENERIC_MSG(WARNING_LEVEL, CERR_STREAM, YELLOW_COLOR, msg)
#define VERBOSE_ERROR(msg)                           VERBOSE_GENERIC_END(ERROR_LEVEL,   CERR_STREAM, RED_COLOR,    msg)
#define VERBOSE_ASSERT(test,msg)  {if (not (test)) { VERBOSE_GENERIC_END(ERROR_LEVEL,   CERR_STREAM, RED_COLOR,    msg) }};





#endif /* ANTMAN_SRC_MIXTURE_CPP_VERBOSE_HPP_ */
