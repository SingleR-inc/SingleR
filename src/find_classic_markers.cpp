#include "Rcpp.h"
#include "utils.h" // must be before raticate, singlepp includes.

#include "singlepp/singlepp.hpp"
#include "raticate/raticate.hpp"

#include <vector>

//[[Rcpp::export(rng=false)]]
Rcpp::List find_classic_markers(int nlabels, int ngenes, Rcpp::List labels, Rcpp::List ref, int de_n, int nthreads) {
    // Storing references to the R objects and creating pointers.
    size_t nref = ref.size();
    if (nref != labels.size()) {
        throw std::runtime_error("'ref' and 'labels' should have the same length");
    }

    std::vector<raticate::Parsed<double, int> > ref_vec;
    ref_vec.reserve(nref);
    std::vector<const tatami::Matrix<double, int>*> ref_ptrs;
    ref_ptrs.reserve(nref);
    std::vector<Rcpp::IntegerVector> lab_vec;
    lab_vec.reserve(nref);
    std::vector<const int*> lab_ptrs;
    lab_ptrs.reserve(nref);

    for (size_t r = 0; r < nref; ++r) {
        ref_vec.push_back(raticate::parse<double, int>(ref[r], true));
        if (ref_vec.back().matrix->nrow() != ngenes) {
            throw std::runtime_error("each entry of 'ref' should have number of rows equal to 'ngenes'");
        }
        ref_ptrs.push_back(ref_vec.back().matrix.get());

        lab_vec.emplace_back(labels[r]);
        if (lab_vec.back().size() != ref_vec.back().matrix->ncol()) {
            throw std::runtime_error("each entry of 'labels' should have length equal to the number of columns in the corresponding entry of 'ref'");
        }
        lab_ptrs.push_back(static_cast<const int*>(lab_vec.back().begin()));
    }
    
    singlepp::ChooseClassicMarkers mrk;
    mrk.set_number(de_n).set_num_threads(nthreads);
    auto store = mrk.run(ref_ptrs, lab_ptrs);

    // Returning everything in R space.
    Rcpp::List output(nlabels);
    for (int l = 0; l < nlabels; ++l) {
        const auto& src = store[l];
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
