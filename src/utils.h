#ifndef UTILS_H 
#define UTILS_H 

#include "Rcpp.h"
#include "singlepp/SinglePP.hpp"
#include "singlepp/IntegratedBuilder.hpp"
#include "tatami/tatami.hpp"
#include "tatami/ext/ArrayView.hpp"
#include <memory>
#include <vector>

typedef Rcpp::XPtr<singlepp::SinglePP::Prebuilt> PrebuiltXPtr;

typedef Rcpp::XPtr<std::vector<singlepp::IntegratedReference> > IntegratedXPtr;

inline std::shared_ptr<tatami::NumericMatrix> tatamize_input(Rcpp::NumericMatrix input) {
    // TODO: allow this to split up depending on the type of input.
    tatami::ArrayView<double> view(static_cast<const double*>(input.begin()), input.size());
    return std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseColumnMatrix<double, int, decltype(view)>(input.nrow(), input.ncol(), std::move(view)));
}

#endif
