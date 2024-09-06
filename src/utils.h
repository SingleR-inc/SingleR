#ifndef UTILS_H 
#define UTILS_H 

#include "Rcpp.h"
#include "Rtatami.h" // before singlepp includes to ensure the tatami_r::parallelize() override is set.
#include "singlepp/singlepp.hpp"

typedef singlepp::TrainedSingleIntersect<int, double> TrainedSingleIntersect;

typedef Rcpp::XPtr<TrainedSingleIntersect> TrainedSingleIntersectPointer;

typedef singlepp::TrainedIntegrated<int> TrainedIntegrated;

typedef Rcpp::XPtr<TrainedIntegrated> TrainedIntegratedPointer;

#endif
