#include "utils.h" // must be before all other includes.

#include <vector>
#include <memory>

//[[Rcpp::export(rng=false)]]
SEXP train_integrated(Rcpp::List test_features, Rcpp::List references, Rcpp::List ref_features, Rcpp::List labels, Rcpp::List prebuilt, int nthreads) {
    size_t nrefs = references.size();

    std::vector<singlepp::TrainIntegratedInput<double, int, int> > inputs;
    inputs.reserve(nrefs);
    std::vector<singlepp::Intersection<int> > intersections(nrefs);
    std::vector<Rcpp::IntegerVector> holding_labs(nrefs);

    for (size_t r = 0; r < nrefs; ++r) {
        Rcpp::RObject curref(references[r]);
        Rtatami::BoundNumericPointer parsed(curref);

        Rcpp::IntegerVector test_ids(test_features[r]);
        Rcpp::IntegerVector ref_ids(ref_features[r]);
        size_t ninter = test_ids.size();
        if (ninter != static_cast<size_t>(ref_ids.size())) {
            throw std::runtime_error("length of each entry of 'test_features' and 'ref_features' should be the same");
        }
        auto& curinter = intersections[r];
        for (size_t i = 0; i < ninter; ++i) {
            curinter.emplace_back(test_ids[i], ref_ids[i]);
        }

        holding_labs[r] = labels[r];
        Rcpp::RObject built = prebuilt[r];
        TrainedSingleIntersectPointer curbuilt(built);

        inputs.push_back(singlepp::prepare_integrated_input_intersect(
            curinter,
            *(parsed->ptr),
            static_cast<const int*>(holding_labs[r].begin()), 
            *curbuilt
        ));
    }

    singlepp::TrainIntegratedOptions opts;
    opts.num_threads = nthreads;
    auto finished = singlepp::train_integrated(std::move(inputs), opts);

    return TrainedIntegratedPointer(new TrainedIntegrated(std::move(finished)), true);
}
