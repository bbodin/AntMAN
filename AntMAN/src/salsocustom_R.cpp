#ifdef HAS_RCPP
#include "salsocustom.h"

// [[Rcpp::export(name= ".salso")]]
Rcpp::List salsoRcpp (const Rcpp::NumericMatrix& eam, int maxClusts,  double Const_Binder, int batchSize, int nScans, int maxThreads, int timeLimit) {
    salso_result_t tmpResult = salsoCpp(Rcpp::as<arma::mat>(eam), maxClusts, Const_Binder, batchSize, nScans, maxThreads, timeLimit);
    return Rcpp::List::create(Rcpp::_["Labels"] = tmpResult.labels, 
                            Rcpp::_["BinderLoss"] = tmpResult.binderLoss,
                            Rcpp::_["NumClusts"] = tmpResult.numClusts,
                            Rcpp::_["NumPermutations"] = tmpResult.nIters, 
                            Rcpp::_["WallClockTime"] = tmpResult.wallClockTime,
                            Rcpp::_["TimeLimitReached"] = tmpResult.timeLimitReached,
                            Rcpp::_["NumThreads"] = tmpResult.numThreads);
}

// [[Rcpp::export(name= ".computeBinderLoss")]]
double computeBinderLossRcpp (const Rcpp::NumericMatrix& eam, const Rcpp::IntegerVector& partitionLabels, double Const_Binder) {
    return computeBinderLossCpp(Rcpp::as<arma::mat>(eam), Rcpp::as<arma::ivec>(partitionLabels), Const_Binder);
}

#endif