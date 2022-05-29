#include "Rcpp.h"
#include "singlepp/IntegratedBuilder.hpp"
#include "utils.h"
#include <vector>
#include <memory>
#include <iostream>

//[[Rcpp::export(rng=false)]]
SEXP integrate_build(Rcpp::IntegerVector test_features, Rcpp::List references, Rcpp::List ref_ids, Rcpp::List labels, Rcpp::List prebuilt) {
    singlepp::IntegratedBuilder builder;

    size_t nrefs = references.size();
    std::vector<std::shared_ptr<tatami::NumericMatrix> > holding_ptrs;
    holding_ptrs.reserve(nrefs);
    std::vector<Rcpp::IntegerVector> holding_labs;
    holding_labs.reserve(nrefs);

    for (size_t r = 0; r < nrefs; ++r) {
        Rcpp::NumericMatrix curref(references[r]);
        holding_ptrs.push_back(tatamize_input(curref));

        Rcpp::IntegerVector curids(ref_ids[r]);
        holding_labs.emplace_back(labels[r]);
        PrebuiltXPtr curbuilt(prebuilt[r]);

        builder.add(
            test_features.size(),
            static_cast<const int*>(test_features.begin()),
            holding_ptrs.back().get(), 
            static_cast<const int*>(curids.begin()),
            static_cast<const int*>(holding_labs.back().begin()), 
            *curbuilt
        );
    }

    auto finished = builder.finish();
    return IntegratedXPtr(new std::vector<singlepp::IntegratedReference>(std::move(finished)), true);
}
