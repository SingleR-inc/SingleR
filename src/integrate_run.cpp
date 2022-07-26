#include "Rcpp.h"

#include "utils.h" // must be before raticate, singlepp includes.

#include "singlepp/singlepp.hpp"
#include "raticate/raticate.hpp"

#include <vector>

//[[Rcpp::export(rng=false)]]
SEXP integrate_run(Rcpp::RObject test, Rcpp::List results, SEXP integrated_build, double quantile, int nthreads) {
    auto curtest = raticate::parse<double, int>(test, true);
    IntegratedXPtr iptr(integrated_build);

    // Setting up the previous results.
    std::vector<Rcpp::IntegerVector> previous_results_vec;
    previous_results_vec.reserve(results.size());
    for (size_t r = 0; r < results.size(); ++r) {
        previous_results_vec.emplace_back(results[r]);
    }

    std::vector<const int*> previous_results;
    previous_results.reserve(results.size());
    for (size_t r = 0; r < results.size(); ++r) {
        previous_results.push_back(static_cast<const int*>(previous_results_vec[r].begin()));
    }

    // Setting up outputs.
    size_t ncells = curtest.matrix->ncol();
    Rcpp::IntegerVector best(ncells);
    Rcpp::NumericVector delta(ncells);

    size_t nrefs = iptr->size();
    Rcpp::NumericMatrix scores(ncells, nrefs);
    std::vector<double*> scores_ptr(nrefs);
    if (nrefs) {
        scores_ptr[0] = static_cast<double*>(scores.begin());
        for (size_t l = 1; l < nrefs; ++l) {
            scores_ptr[l] = scores_ptr[l-1] + ncells;
        }
    }

    // Running the integrated scoring.
    singlepp::IntegratedScorer scorer;
    scorer.set_num_threads(nthreads).set_quantile(quantile);
    scorer.run(
        curtest.matrix.get(), 
        previous_results, 
        *iptr, 
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

