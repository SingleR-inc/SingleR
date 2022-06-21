#include "Rcpp.h"
#include <vector>
#include "utils.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List find_classic_markers(int nlabels, int ngenes, Rcpp::IntegerVector left, Rcpp::IntegerVector right, Rcpp::List ref, Rcpp::List available, int de_n, int nthreads) {
    // Storing references to the R objects and creating pointers.
    size_t nref = ref.size();
    if (nref != available.size()) {
        throw std::runtime_error("'ref' and 'available' should have the same length");
    }

    std::vector<Rcpp::NumericMatrix> ref_vec;
    ref_vec.reserve(nref);
    std::vector<Rcpp::LogicalVector> avail_vec;
    avail_vec.reserve(nref);
    for (size_t r = 0; r < nref; ++r) {
        ref_vec.emplace_back(ref[r]);
        if (ref_vec.back().rows() != ngenes) {
            throw std::runtime_error("each entry of 'ref' should have number of rows equal to 'ngenes'");
        }

        avail_vec.emplace_back(available[r]);
        if (avail_vec.back().size() != ref_vec.back().cols()) {
            throw std::runtime_error("each entry of 'available' should have length equal to the number of columns in the corresponding entry of 'ref'");
        }
    }
    
    auto left_ptr = static_cast<const int*>(left.begin());
    auto right_ptr = static_cast<const int*>(right.begin());
    size_t npairs = left.size();
    if (right.size() != npairs) {
        throw std::runtime_error("'left' and 'right' should have the same length");
    }
    for (size_t p = 0; p < npairs; ++p) {
        if (left_ptr[p] < 0 || left_ptr[p] >= nlabels) {
            throw std::runtime_error("label indexes of 'left' are out of range");
        }
        if (right_ptr[p] < 0 || right_ptr[p] >= nlabels) {
            throw std::runtime_error("label indexes of 'right' are out of range");
        }
    }

    std::vector<const double*> ref_ptrs;
    ref_ptrs.reserve(nref);
    std::vector<const int*> avail_ptrs;
    avail_ptrs.reserve(nref);
    for (size_t r = 0; r < nref; ++r) {
        ref_ptrs.push_back(ref_vec[r].begin());
        avail_ptrs.push_back(avail_vec[r].begin());
    }

    if (de_n > ngenes) {
        de_n = ngenes;
    }

    // Setting up the stores.
    std::vector<std::vector<std::vector<int> > > store(nlabels, std::vector<std::vector<int> >(nlabels));

    num_threads = nthreads;
    parallelize(npairs, [&](size_t first, size_t last) -> void {
        std::vector<std::pair<double, int> > sorter(ngenes);
        for (size_t p = first; p < last; ++p) {
            auto curleft = left_ptr[p];
            auto curright = right_ptr[p];

            auto sIt = sorter.begin();
            for (int g = 0; g < ngenes; ++g, ++sIt) {
                sIt->first = 0;
                sIt->second = g;
            }

            for (size_t r = 0; r < nref; ++r) {
                auto cur_avail = avail_ptrs[r];
                if (!cur_avail[curleft] || !cur_avail[curright]) {
                    continue;
                }

                auto leftcol = ref_ptrs[r] + ngenes * curleft;
                auto rightcol = ref_ptrs[r] + ngenes * curright;
                auto sIt = sorter.begin();
                for (int g = 0; g < ngenes; ++g, ++leftcol, ++rightcol, ++sIt) {
                    sIt->first += -(*leftcol - *rightcol); // reversing log-fold change, so we sort for the largest logFCs.
                }
            }

            // partial sort is guaranteed to be stable due to the second index resolving ties.
            std::partial_sort(sorter.begin(), sorter.begin() + de_n, sorter.end());

            std::vector<int> stuff;
            stuff.reserve(de_n);
            for (int g = 0; g < de_n && sorter[g].first < 0; ++g) { // only keeping those with positive log-fold changes (negative after reversing).
                stuff.push_back(sorter[g].second + 1); // for 1-based indexing.
            }
            store[curleft][curright] = std::move(stuff);
        }
    });

    // Returning everything in R space.
    Rcpp::List output(nlabels);
    for (int l = 0; l < nlabels; ++l) {
        const auto& src = store[l];
        Rcpp::List dest(nlabels);
        for (int l2 = 0; l2 < nlabels; ++l2) {
            dest[l2] = Rcpp::IntegerVector(src[l2].begin(), src[l2].end());
        }
        output[l] = dest;
    }

    return output;
}
