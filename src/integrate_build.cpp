#include "Rcpp.h"
#include "singlepp/IntegratedBuilder.hpp"
#include "utils.h"
#include <vector>
#include <memory>

//[[Rcpp::export(rng=false)]]
SEXP integrate_build(Rcpp::IntegerVector test_features, Rcpp::List references, Rcpp::List ref_ids, Rcpp::List labels, Rcpp::List prebuilt) {
    singlepp::IntegratedBuilder builder;

    size_t nrefs = references.size();
    for (size_t r = 0; r < nrefs; ++r) {
        Rcpp::NumericMatrix curref(references[r]);
        auto curptr = tatamize_input(curref);
        Rcpp::IntegerVector curids(ref_ids[r]);
        Rcpp::IntegerVector curlab(labels[r]);
        PrebuiltXPtr curbuilt(prebuilt[r]);

        builder.add(
            test_features.size(),
            static_cast<const int*>(test_features.begin()),
            curptr.get(), 
            static_cast<const int*>(curids.begin()),
            static_cast<const int*>(curlab.begin()), 
            *curbuilt
        );
    }

    return IntegratedXPtr(new std::vector<singlepp::IntegratedReference>(std::move(builder.finish())), true);
}
