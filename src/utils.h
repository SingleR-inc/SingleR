#ifndef UTILS_H 
#define UTILS_H 

#include "Rcpp.h"
#include "singlepp/SinglePP.hpp"
#include "singlepp/IntegratedBuilder.hpp"
#include <memory>
#include <vector>

typedef Rcpp::XPtr<singlepp::SinglePP::Prebuilt> PrebuiltXPtr;

typedef Rcpp::XPtr<std::vector<singlepp::IntegratedReference> > IntegratedXPtr;

#endif
