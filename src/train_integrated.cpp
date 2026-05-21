#include "utils.h" // must be before all other includes.

#include <vector>
#include <memory>
#include <cstddef>

//[[Rcpp::export(rng=false)]]
SEXP train_integrated(
    int test_nrow,
    Rcpp::List test_features,
    Rcpp::List references,
    Rcpp::List ref_features,
    Rcpp::List labels,
    Rcpp::List markers,
    int nthreads
) {
    const std::size_t nrefs = references.size();
    std::vector<singlepp::TrainIntegratedInput<double, int, int> > inputs;
    inputs.reserve(nrefs);
    std::vector<Rcpp::IntegerVector> holding_labs(nrefs);

    for (std::size_t r = 0; r < nrefs; ++r) {
        Rcpp::RObject curref(references[r]);
        Rtatami::BoundNumericPointer parsed(curref);

        Rcpp::IntegerVector test_ids(test_features[r]);
        Rcpp::IntegerVector ref_ids(ref_features[r]);
        const std::size_t ninter = test_ids.size();
        if (ninter != static_cast<std::size_t>(ref_ids.size())) {
            throw std::runtime_error("length of each entry of 'test_features' and 'ref_features' should be the same");
        }

        singlepp::Intersection<int> curinter;
        curinter.reserve(ninter);
        for (std::size_t i = 0; i < ninter; ++i) {
            curinter.emplace_back(test_ids[i], ref_ids[i]);
        }

        // Setting up the markers. We assume that these are already 0-indexed on the R side.
        Rcpp::List my_markers = markers[r];
        const std::size_t ngroups = my_markers.size();
        singlepp::PerLabelMarkers<int> converted_markers(ngroups);
        for (std::size_t m = 0; m < ngroups; ++m) {
            Rcpp::IntegerVector curmarkers(my_markers[m]);
            auto& dest = converted_markers[m];
            dest.insert(dest.end(), curmarkers.begin(), curmarkers.end());
        }

        holding_labs[r] = labels[r]; // holding a reference to avoid GC of the array, if it was realized from an ALTREP.
        inputs.push_back(singlepp::prepare_integrated_input<int, double, int>(
            test_nrow,
            std::move(curinter),
            parsed->ptr,
            static_cast<const int*>(holding_labs[r].begin()), 
            std::move(converted_markers)
        ));
    }

    singlepp::TrainIntegratedOptions opts;
    opts.num_threads = nthreads;
    auto finished = singlepp::train_integrated(std::move(inputs), opts);

    return TrainedIntegratedPointer(new TrainedIntegrated(std::move(finished)), true);
}
