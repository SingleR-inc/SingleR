#include "utils.h" // must be before all other includes.

#include <vector>
#include <memory>

//[[Rcpp::export(rng=false)]]
SEXP train_integrated(Rcpp::IntegerVector test_features, Rcpp::List references, Rcpp::List ref_ids, Rcpp::List labels, Rcpp::List prebuilt, int nthreads) {
    size_t nrefs = references.size();

    std::vector<singlepp::TrainIntegratedInput<double, int, int> > inputs;
    inputs.reserve(nrefs);
    std::vector<Rcpp::IntegerVector> holding_labs;
    holding_labs.reserve(nrefs);

    for (size_t r = 0; r < nrefs; ++r) {
        Rcpp::RObject curref(references[r]);
        Rtatami::BoundNumericPointer parsed(curref);

        Rcpp::IntegerVector curids(ref_ids[r]);
        holding_labs.emplace_back(labels[r]);
        Rcpp::RObject built = prebuilt[r];
        TrainedSingleIntersectPointer curbuilt(built);

        inputs.push_back(singlepp::prepare_integrated_input_intersect(
            static_cast<int>(test_features.size()),
            static_cast<const int*>(test_features.begin()),
            *(parsed->ptr),
            static_cast<const int*>(curids.begin()),
            static_cast<const int*>(holding_labs.back().begin()), 
            *curbuilt
        ));
    }

    singlepp::TrainIntegratedOptions opts;
    opts.num_threads = nthreads;
    auto finished = singlepp::train_integrated(std::move(inputs), opts);

    return TrainedIntegratedPointer(new TrainedIntegrated(std::move(finished)), true);
}
