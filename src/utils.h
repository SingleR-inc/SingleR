#ifndef UTILS_H 
#define UTILS_H 

#include "Rcpp.h"
#include "Rtatami.h"

// must be before singlepp includes.
#define SINGLEPP_CUSTOM_PARALLEL tatami_r::parallelize 
#define __ERROR_PRINTER_OVERRIDE__  REprintf // avoid R CMD check warnings about stderr in Annoy.

#include "singlepp/singlepp.hpp"

typedef Rcpp::XPtr<singlepp::BasicBuilder::Prebuilt> PrebuiltXPtr;

typedef Rcpp::XPtr<singlepp::IntegratedReferences> IntegratedXPtr;

#endif
