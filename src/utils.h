#ifndef UTILS_H 
#define UTILS_H 

#include <cmath>
#include <vector>
#include <thread>
#include <memory>
#include "Rcpp.h"

extern int num_threads;

#define RATICATE_PARALLELIZE_UNKNOWN // must be before raticate includes.
#include "raticate/raticate.hpp"

template<class Function>
void parallelize(size_t n, Function function, int num_threads) {
    raticate::parallelize<double, int>(n, function, num_threads);
}

// must be before singlepp includes.
#define SINGLEPP_CUSTOM_PARALLEL parallelize 
#define __ERROR_PRINTER_OVERRIDE__  REprintf // avoid R CMD check warnings about stderr in Annoy.

#include "singlepp/singlepp.hpp"

typedef Rcpp::XPtr<singlepp::BasicBuilder::Prebuilt> PrebuiltXPtr;

typedef Rcpp::XPtr<std::vector<singlepp::IntegratedReference> > IntegratedXPtr;

#endif
