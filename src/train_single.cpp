#include "utils.h" 
#include "BiocNeighbors.h"

#include <vector>
#include <memory>

//' @importFrom Rcpp sourceCpp
//' @useDynLib SingleR
//[[Rcpp::export(rng=false)]]
SEXP train_single(Rcpp::IntegerVector test_features, Rcpp::RObject ref, Rcpp::IntegerVector ref_features, Rcpp::IntegerVector labels, Rcpp::List markers, Rcpp::RObject builder, int nthreads) {
    singlepp::TrainSingleOptions opts;
    opts.num_threads = nthreads;
    opts.top = -1; // Use all available markers; assume subsetting was applied on the R side.

    BiocNeighbors::BuilderPointer bptr(builder);
    opts.trainer = std::shared_ptr<BiocNeighbors::Builder>(std::shared_ptr<BiocNeighbors::Builder>{}, bptr.get()); // make a no-op shared pointer.

    Rtatami::BoundNumericPointer parsed(ref);
    int NR = parsed->ptr->nrow();
    int NC = parsed->ptr->ncol();
    if (static_cast<int>(labels.size()) != NC) {
        throw std::runtime_error("length of 'labels' is equal to the number of columns of 'ref'");
    }

    // Setting up the markers. We assume that these are already 0-indexed on the R side.
    size_t ngroups = markers.size();
    singlepp::Markers<int> markers2(ngroups);
    for (size_t m = 0; m < ngroups; ++m) {
        Rcpp::List curmarkers(markers[m]);
        auto& curmarkers2 = markers2[m];
        size_t inner_ngroups = curmarkers.size();
        curmarkers2.resize(inner_ngroups);

        for (size_t n = 0; n < inner_ngroups; ++n) {
            Rcpp::IntegerVector seq(curmarkers[n]);
            auto& seq2 = curmarkers2[n];
            seq2.insert(seq2.end(), seq.begin(), seq.end());
        }
    }

    // Preparing the features.
    size_t ninter = test_features.size();
    if (ninter != static_cast<size_t>(ref_features.size())) {
        throw std::runtime_error("length of 'test_features' and 'ref_features' should be the same");
    }
    singlepp::Intersection<int> inter;
    inter.reserve(test_features.size());
    for (size_t i = 0; i < ninter; ++i) {
        inter.emplace_back(test_features[i], ref_features[i]);
    }

    // Building the indices.
    auto built = singlepp::train_single_intersect(
        inter,
        *(parsed->ptr),
        static_cast<const int*>(labels.begin()),
        std::move(markers2),
        opts
    );

    return TrainedSingleIntersectPointer(new TrainedSingleIntersect(std::move(built)), true);
}

//[[Rcpp::export(rng=false)]]
Rcpp::IntegerVector get_ref_subset(SEXP built) {
    TrainedSingleIntersectPointer ptr(built);
    const auto& rsub = ptr->get_ref_subset();
    return Rcpp::IntegerVector(rsub.begin(), rsub.end());
}

//[[Rcpp::export(rng=false)]]
Rcpp::LogicalVector is_valid_built(SEXP built) {
    return Rf_ScalarLogical(!!R_ExternalPtrAddr(built));
}
