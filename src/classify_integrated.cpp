#include "utils.h" // must be before all other includes.

#include <vector>

//[[Rcpp::export(rng=false)]]
SEXP classify_integrated(Rcpp::RObject test, Rcpp::List results, SEXP integrated_build, double quantile, int nthreads) {
    Rtatami::BoundNumericPointer curtest(test);
    TrainedIntegratedPointer iptr(integrated_build);

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
    size_t ncells = curtest->ptr->ncol();
    Rcpp::IntegerVector best(ncells);
    Rcpp::NumericVector delta(ncells);

    singlepp::ClassifyIntegratedBuffers<int, double> buffers;
    buffers.best = static_cast<int*>(best.begin());
    buffers.delta = static_cast<double*>(delta.begin());

    size_t nrefs = iptr->num_references();
    Rcpp::NumericMatrix scores(ncells, nrefs);
    if (nrefs) {
        buffers.scores.resize(nrefs);
        buffers.scores[0] = static_cast<double*>(scores.begin());
        for (size_t l = 1; l < nrefs; ++l) {
            buffers.scores[l] = buffers.scores[l-1] + ncells;
        }
    }

    // Running the integrated scoring.
    singlepp::ClassifyIntegratedOptions<double> opts;
    opts.num_threads = nthreads;
    opts.quantile = quantile;
    singlepp::classify_integrated(*(curtest->ptr), previous_results, *iptr, buffers, opts);

    return Rcpp::List::create(
        Rcpp::Named("best") = best,
        Rcpp::Named("scores") = scores,
        Rcpp::Named("delta") = delta
    );
}

