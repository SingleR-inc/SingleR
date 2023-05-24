#include "utils.h" // must be before all other includes.

#include <vector>
#include <memory>

//[[Rcpp::export(rng=false)]]
SEXP integrate_build(Rcpp::IntegerVector test_features, Rcpp::List references, Rcpp::List ref_ids, Rcpp::List labels, Rcpp::List prebuilt, int nthreads) {
    singlepp::IntegratedBuilder builder;
    builder.set_num_threads(nthreads);

    size_t nrefs = references.size();
    std::vector<Rcpp::IntegerVector> holding_labs;
    holding_labs.reserve(nrefs);

    for (size_t r = 0; r < nrefs; ++r) {
        Rcpp::RObject curref(references[r]);
        auto parsed = Rtatami::BoundNumericPointer(curref);

        Rcpp::IntegerVector curids(ref_ids[r]);
        holding_labs.emplace_back(labels[r]);
        Rcpp::RObject built = prebuilt[r];
        PrebuiltXPtr curbuilt(built);

        builder.add(
            test_features.size(),
            static_cast<const int*>(test_features.begin()),
            parsed->ptr.get(),
            static_cast<const int*>(curids.begin()),
            static_cast<const int*>(holding_labs.back().begin()), 
            *curbuilt
        );
    }

    auto finished = builder.finish();
    return IntegratedXPtr(new singlepp::IntegratedReferences(std::move(finished)), true);
}
