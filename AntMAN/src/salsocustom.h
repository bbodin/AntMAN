#pragma once
#include <algorithm>
#include <numeric>
#include <vector>
#include <random>
#include <limits>
#include <chrono>

#include "utils.h"
#ifndef HAS_RCPP
#include <armadillo>
#define ARMA_NO_DEBUG
#else
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

typedef arma::uword ind_t; // type for all counters, indices, and labels
typedef arma::Col<ind_t> labelVec_t; // label vectors
static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 required");
extern double negative_infinity;

struct salso_result_t {
  ind_t nIters;
  unsigned int wallClockTime;
  bool timeLimitReached;
  int numThreads;
  labelVec_t labels;
	ind_t numClusts;
	double binderLoss;
	salso_result_t (ind_t numElems) : nIters(0), wallClockTime(0), timeLimitReached(false), numThreads(0), labels (numElems), numClusts (0), binderLoss (negative_infinity) {}
};
salso_result_t salsoCpp (const arma::mat& eam, ind_t maxClusts, double Const_Binder, ind_t batchSize, ind_t nScans, unsigned int maxThreads, unsigned int timeLimit);

double computeBinderLossCpp (const arma::mat& eam, const arma::ivec& partitionLabels, double Const_Binder);
std::vector<ind_t> randperm (ind_t N);
