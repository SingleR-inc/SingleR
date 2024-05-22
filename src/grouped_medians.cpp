#include "utils.h" // must be before all other includes.

#include "tatami_stats/tatami_stats.hpp"
#include <vector>

//[[Rcpp::export(rng=false)]]
Rcpp::NumericMatrix grouped_medians(Rcpp::RObject ref, Rcpp::IntegerVector groups, int ngroups, int nthreads) {
    Rtatami::BoundNumericPointer parsed(ref);
    size_t NR = parsed->ptr->nrow();

    Rcpp::NumericMatrix output(NR, ngroups);
    double* optr = static_cast<double*>(output.begin());
    std::vector<double*> optrs;
    optrs.reserve(ngroups);
    for (int x = 0; x < ngroups; ++x, optr += NR) {
        optrs.emplace_back(optr);
    }

    std::vector<int> group_sizes(ngroups);
    for (auto g : groups) {
        ++(group_sizes[g]);
    }

    tatami_stats::grouped_medians::Options opt;
    opt.num_threads = nthreads;
    tatami_stats::grouped_medians::apply(
        /* byrow = */ true,
        parsed->ptr.get(), 
        static_cast<const int*>(groups.begin()), 
        group_sizes,
        optrs.data(),
        opt
    );
    return output;
}
