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

    // Setting up the markers. We assume that these are already 0-indexed on the R side.
    singlepp::Markers<int> markers2(markers.size());
    for (size_t m = 0; m < markers.size(); ++m) {
        Rcpp::List curmarkers(markers[m]);
        auto& curmarkers2 = markers2[m];
        curmarkers2.resize(curmarkers.size());

        for (size_t n = 0; n < curmarkers.size(); ++n) {
            Rcpp::IntegerVector seq(curmarkers[n]);
            auto& seq2 =  curmarkers2[n];
            seq2.insert(seq2.end(), seq.begin(), seq.end());
        }
    }

    // Building the indices.
    Rtatami::BoundNumericPointer parsed(ref);
    auto built = singlepp::train_single_intersect(
        static_cast<int>(test_features.size()),
        static_cast<const int*>(test_features.begin()),
        *(parsed->ptr),
        static_cast<const int*>(ref_features.begin()),
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
