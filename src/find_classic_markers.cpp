#include "utils.h" // must be before all other includes.

#include <vector>

//[[Rcpp::export(rng=false)]]
Rcpp::List find_classic_markers(int nlabels, int ngenes, Rcpp::List labels, Rcpp::List ref, int de_n, int nthreads) {
    // Storing references to the R objects and creating pointers.
    size_t nref = ref.size();
    if (nref != static_cast<size_t>(labels.size())) {
        throw std::runtime_error("'ref' and 'labels' should have the same length");
    }

    std::vector<const tatami::Matrix<double, int>*> ref_ptrs;
    ref_ptrs.reserve(nref);
    std::vector<Rcpp::IntegerVector> lab_vec;
    lab_vec.reserve(nref);
    std::vector<const int*> lab_ptrs;
    lab_ptrs.reserve(nref);

    for (size_t r = 0; r < nref; ++r) {
        Rcpp::RObject current = ref[r];
        Rtatami::BoundNumericPointer parsed(current);
        const auto& ptr = parsed->ptr;

        if (ptr->nrow() != ngenes) {
            throw std::runtime_error("each entry of 'ref' should have number of rows equal to 'ngenes'");
        }
        ref_ptrs.push_back(ptr.get());

        lab_vec.emplace_back(labels[r]);
        if (lab_vec.back().size() != ptr->ncol()) {
            throw std::runtime_error("each entry of 'labels' should have length equal to the number of columns in the corresponding entry of 'ref'");
        }
        lab_ptrs.push_back(static_cast<const int*>(lab_vec.back().begin()));
    }
    
    singlepp::ChooseClassicMarkersOptions opts;
    opts.number = de_n;
    opts.num_threads = nthreads;
    auto store = singlepp::choose_classic_markers(ref_ptrs, lab_ptrs, opts);

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
