#include "Rcpp.h"

#include "utils.h" // must be before raticate, singlepp includes.

#include "singlepp/singlepp.hpp"
#include "raticate/raticate.hpp"

#include <vector>

//' @importFrom Rcpp sourceCpp
//' @useDynLib SingleR
//[[Rcpp::export(rng=false)]]
SEXP run(Rcpp::RObject test, Rcpp::IntegerVector subset, SEXP prebuilt, double quantile, bool use_fine_tune, double fine_tune_threshold, int nthreads) {
    num_threads = nthreads;
    auto parsed = raticate::parse<double, int>(test, true);
    PrebuiltXPtr built(prebuilt);

    // Setting up outputs.
    size_t ncells = parsed.matrix->ncol();
    Rcpp::IntegerVector best(ncells);
    Rcpp::NumericVector delta(ncells);

    size_t nlabels = built->num_labels();
    Rcpp::NumericMatrix scores(ncells, nlabels);
    std::vector<double*> scores_ptr(nlabels);
    if (nlabels) {
        scores_ptr[0] = static_cast<double*>(scores.begin());
        for (size_t l = 1; l < nlabels; ++l) {
            scores_ptr[l] = scores_ptr[l-1] + ncells;
        }
    }

    // Running the analysis.
    singlepp::Classifier runner;
    runner.set_quantile(quantile).set_fine_tune(use_fine_tune).set_fine_tune_threshold(fine_tune_threshold);
    runner.run(
        parsed.matrix.get(), 
        *built, 
        static_cast<const int*>(subset.begin()),
        static_cast<int*>(best.begin()),
        scores_ptr,
        static_cast<double*>(delta.begin())
    );

    return Rcpp::List::create(
        Rcpp::Named("best") = best,
        Rcpp::Named("scores") = scores,
        Rcpp::Named("delta") = delta
    );
}
