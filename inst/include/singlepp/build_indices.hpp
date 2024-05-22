#ifndef SINGLEPP_BUILD_INDICES_HPP
#define SINGLEPP_BUILD_INDICES_HPP

#include "macros.hpp"

#include "knncolle/knncolle.hpp"
#include "tatami/tatami.hpp"

#include "process_features.hpp"
#include "scaled_ranks.hpp"

#include <vector>
#include <memory>
#include <algorithm>

namespace singlepp {

inline size_t get_nlabels(size_t n, const int* labels) { 
    if (n == 0) {
        throw std::runtime_error("reference dataset must have at least one column");
    }
    return *std::max_element(labels, labels + n) + 1;
}

struct Reference {
    std::vector<RankedVector<int, int> > ranked;
    std::shared_ptr<knncolle::Base<int, double> > index;
};

template<class Builder>
std::vector<Reference> build_indices(const tatami::Matrix<double, int>* ref, const int* labels, const std::vector<int>& subset, const Builder& build, int nthreads) {
    size_t NC = ref->ncol();
    size_t nlabels = get_nlabels(NC, labels);
    std::vector<int> label_count(nlabels);
    for (size_t i = 0; i < NC; ++i) {
        ++label_count[labels[i]];
    }

    size_t NR = subset.size();
    std::vector<Reference> nnrefs(nlabels);
    std::vector<std::vector<double> > nndata(nlabels);
    for (size_t l = 0; l < nlabels; ++l) {
        if (label_count[l] == 0) {
            throw std::runtime_error(std::string("no entries for label ") + std::to_string(l));
        }
        nnrefs[l].ranked.resize(label_count[l]);
        nndata[l].resize(label_count[l] * NR);
    }

    std::vector<size_t> offsets(NC);
    {
        std::vector<size_t> counter(nlabels);
        for (size_t c = 0; c < NC; ++c) {
            auto& o = counter[labels[c]];
            offsets[c] = o;
            ++o;
        }
    }

    SubsetSorter subsorter(subset);

    tatami::parallelize([&](int, int start, int len) -> void {
        RankedVector<double, int> ranked(NR);
        std::vector<double> buffer(ref->nrow());
        auto wrk = tatami::consecutive_extractor<false>(ref, false, start, len, subsorter.extraction_subset());

        for (int c = start, end = start + len; c < end; ++c) {
            auto ptr = wrk->fetch(c, buffer.data());
            subsorter.fill_ranks(ptr, ranked); 

            auto curlab = labels[c];
            auto curoff = offsets[c];
            auto scaled = nndata[curlab].data() + curoff * NR;
            scaled_ranks(ranked, scaled); 

            // Storing as a pair of ints to save space; as long
            // as we respect ties, everything should be fine.
            auto& stored_ranks = nnrefs[curlab].ranked[curoff];
            stored_ranks.reserve(ranked.size());
            simplify_ranks(ranked, stored_ranks);
        }
    }, ref->ncol(), nthreads);

#ifndef SINGLEPP_CUSTOM_PARALLEL
    #pragma omp parallel for num_threads(nthreads)
    for (size_t l = 0; l < nlabels; ++l) {
#else
    SINGLEPP_CUSTOM_PARALLEL([&](int, size_t start, size_t len) -> void {
    for (size_t l = start, end = start + len; l < end; ++l) {
#endif
        nnrefs[l].index = build(NR, label_count[l], nndata[l].data());

        // Trying to free the memory as we go along, to offset the copying
        // of nndata into the memory store owned by the knncolle index.
        nndata[l].clear();
        nndata[l].shrink_to_fit();

#ifndef SINGLEPP_CUSTOM_PARALLEL
    }
#else
    }
    }, nlabels, nthreads);
#endif

    return nnrefs;
}

}

#endif
