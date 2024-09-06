#include "utils.h" // must be before raticate, singlepp includes.

#include <vector>

//' @importFrom Rcpp sourceCpp
//' @useDynLib SingleR
//[[Rcpp::export(rng=false)]]
SEXP classify_single(Rcpp::RObject test, SEXP prebuilt, double quantile, bool use_fine_tune, double fine_tune_threshold, int nthreads) {
    Rtatami::BoundNumericPointer parsed(test);
    TrainedSingleIntersectPointer built(prebuilt);

    // Setting up outputs.
    size_t ncells = parsed->ptr->ncol();
    Rcpp::IntegerVector best(ncells);
    Rcpp::NumericVector delta(ncells);

    singlepp::ClassifySingleBuffers buffers;
    buffers.best = static_cast<int*>(best.begin());
    buffers.delta = static_cast<double*>(delta.begin());

    size_t nlabels = built->num_labels();
    Rcpp::NumericMatrix scores(ncells, nlabels);
    if (nlabels) {
        buffers.scores.resize(nlabels);
        buffers.scores[0] = static_cast<double*>(scores.begin());
        for (size_t l = 1; l < nlabels; ++l) {
            buffers.scores[l] = buffers.scores[l-1] + ncells;
        }
    }

    // Running the analysis.
    singlepp::ClassifySingleOptions opts;
    opts.num_threads = nthreads;
    opts.quantile = quantile;
    opts.fine_tune = use_fine_tune;
    opts.fine_tune_threshold = fine_tune_threshold;
    singlepp::classify_single_intersect(*(parsed->ptr), *built, buffers, opts);

    return Rcpp::List::create(
        Rcpp::Named("best") = best,
        Rcpp::Named("scores") = scores,
        Rcpp::Named("delta") = delta
    );
}
