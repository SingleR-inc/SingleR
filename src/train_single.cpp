#include "utils.h" 

#include <vector>
#include <memory>

//' @importFrom Rcpp sourceCpp
//' @useDynLib SingleR
//[[Rcpp::export(rng=false)]]
SEXP train_single(
    int test_nrow, 
    Rcpp::IntegerVector test_features,
    Rcpp::RObject ref,
    Rcpp::IntegerVector ref_features,
    Rcpp::IntegerVector labels,
    Rcpp::List markers,
    int nthreads
) {
    // We use all available markers; assume subsetting was applied on the R side.
    singlepp::TrainSingleOptions opts;
    opts.num_threads = nthreads;

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
    std::vector<int> ref_subset;
    auto built = singlepp::train_single(
        test_nrow,
        inter,
        *(parsed->ptr),
        static_cast<const int*>(labels.begin()),
        std::move(markers2),
        &ref_subset,
        opts
    );

    return Rcpp::List::create(
        Rcpp::Named("index") = TrainedSinglePointer(new TrainedSingle(std::move(built)), true),
        Rcpp::Named("ref_subset") = Rcpp::IntegerVector(ref_subset.begin(), ref_subset.end())
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::LogicalVector is_valid_built(SEXP built) {
    return Rf_ScalarLogical(!!R_ExternalPtrAddr(built));
}
