#include "Rcpp.h"

#include "utils.h" // must be before raticate, singlepp includes.

#include "singlepp/singlepp.hpp"
#include "raticate/raticate.hpp"

#include <vector>
#include <memory>

//[[Rcpp::export(rng=false)]]
SEXP integrate_build(Rcpp::IntegerVector test_features, Rcpp::List references, Rcpp::List ref_ids, Rcpp::List labels, Rcpp::List prebuilt, int nthreads) {
    num_threads = nthreads;
    singlepp::IntegratedBuilder builder;

    size_t nrefs = references.size();
    std::vector<raticate::Parsed<double, int> > holding_mats;
    holding_mats.reserve(nrefs);
    std::vector<Rcpp::IntegerVector> holding_labs;
    holding_labs.reserve(nrefs);

    for (size_t r = 0; r < nrefs; ++r) {
        Rcpp::RObject curref(references[r]);
        auto parsed = raticate::parse<double, int>(curref, true);
        holding_mats.push_back(std::move(parsed));

        Rcpp::IntegerVector curids(ref_ids[r]);
        holding_labs.emplace_back(labels[r]);
        PrebuiltXPtr curbuilt(prebuilt[r]);

        builder.add(
            test_features.size(),
            static_cast<const int*>(test_features.begin()),
            holding_mats.back().matrix.get(), 
            static_cast<const int*>(curids.begin()),
            static_cast<const int*>(holding_labs.back().begin()), 
            *curbuilt
        );
    }

    auto finished = builder.finish();
    return IntegratedXPtr(new std::vector<singlepp::IntegratedReference>(std::move(finished)), true);
}
