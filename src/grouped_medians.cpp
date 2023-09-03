#include "utils.h" // must be before all other includes.

#include <vector>

//[[Rcpp::export(rng=false)]]
Rcpp::NumericMatrix grouped_medians(Rcpp::RObject ref, Rcpp::IntegerVector groups, int ngroups, int nthreads) {
    Rtatami::BoundNumericPointer parsed(ref);
    Rcpp::NumericMatrix output(ngroups, parsed->ptr->nrow());

    std::vector<int> group_sizes(ngroups);
    for (auto g : groups) {
        ++(group_sizes[g]);
    }

    tatami::row_medians_by_group(
        parsed->ptr.get(), 
        static_cast<const int*>(groups.begin()), 
        group_sizes,
        static_cast<double*>(output.begin()), 
        nthreads
    );
    return output;
}
