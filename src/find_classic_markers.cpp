#include "utils.h" // must be before all other includes.

#include "Rcpp.h"

#include "sanisizer/sanisizer.hpp"
#include "singler_classic_markers/singler_classic_markers.hpp"

#include <stdexcept>

//[[Rcpp::export(rng=false)]]
Rcpp::List find_classic_markers(
    Rcpp::RObject mat,
    int nlabels,
    Rcpp::IntegerVector labels,
    Rcpp::Nullable<Rcpp::IntegerVector> blocks,
    Rcpp::Nullable<int> de_n,
    int nthreads
) {
    Rtatami::BoundNumericPointer parsed(mat);
    const auto& matrix = *(parsed->ptr);
    if (!sanisizer::is_equal(matrix.ncol(), labels.size())) {
        throw std::runtime_error("number of columns in 'ref' should equal the length of 'labels'");
    }

    std::optional<std::size_t> number;
    if (!de_n.isNull()) {
        Rcpp::IntegerVector val(de_n);
        if (val.size() != 1 || val[0] < 0) {
            throw std::runtime_error("'de_n' should be a positive scalar if provided");
        }
        number = val[0];
    }

    std::vector<std::vector<std::vector<int> > > markers;

    if (blocks.isNull()) {
        singler_classic_markers::ChooseOptions opts;
        opts.num_threads = nthreads;
        opts.number = number;
        markers = singler_classic_markers::choose_index(matrix, static_cast<const int*>(labels.begin()), opts);

    } else {
        Rcpp::IntegerVector bb(blocks);
        if (!sanisizer::is_equal(matrix.ncol(), bb.size())) {
            throw std::runtime_error("number of columns in 'ref' should equal the length of 'blocks' if provided");
        }
        singler_classic_markers::ChooseBlockedOptions opts;
        opts.num_threads = nthreads;
        opts.number = number;
        markers = singler_classic_markers::choose_blocked_index(matrix, static_cast<const int*>(labels.begin()), static_cast<const int*>(bb.begin()), opts);
    }

    // Returning everything in R space.
    Rcpp::List output(nlabels);
    for (int l = 0; l < nlabels; ++l) {
        const auto& src = markers[l];
        Rcpp::List dest(nlabels);
        for (int l2 = 0; l2 < nlabels; ++l2) {
            Rcpp::IntegerVector current(src[l2].begin(), src[l2].end());
            for (auto& c : current) {
                ++c;
            }
            dest[l2] = std::move(current);
        }
        output[l] = dest;
    }

    return output;
}
