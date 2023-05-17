#ifndef SINGLEPP_SCALED_RANKS_HPP
#define SINGLEPP_SCALED_RANKS_HPP

#include "macros.hpp"

#include <algorithm>
#include <vector>
#include <cmath>
#include <unordered_map>

namespace singlepp {

template<typename Stat, typename Index>
using RankedVector = std::vector<std::pair<Stat, Index> >;

// This class sanitizes any user-provided subsets so that we can provide a
// sorted and unique subset to the tatami extractor. We then undo the sorting
// to use the original indices in the rank filler. This entire thing is
// necessary as the behavior of the subsets isn't something that the user can
// easily control (e.g., if the reference/test datasets do not use the same
// feature ordering, in which case the subset is necessarily unsorted).
struct SubsetSorter {
    bool use_sorted_subset = false;
    const std::vector<int>* original_subset;
    std::vector<int> sorted_subset, original_indices;

    SubsetSorter(const std::vector<int>& sub) : original_subset(&sub) {
        size_t num_subset = sub.size();
        for (size_t i = 1; i < num_subset; ++i) {
            if (sub[i] <= sub[i-1]) {
                use_sorted_subset = true;
                break;
            }
        }

        if (use_sorted_subset) {
            std::vector<std::pair<int, int> > store;
            store.reserve(num_subset);
            for (size_t i = 0; i < num_subset; ++i) {
                store.emplace_back(sub[i], i);
            }
            
            std::sort(store.begin(), store.end());
            sorted_subset.reserve(num_subset);
            original_indices.resize(num_subset);
            for (const auto& s : store) {
                if (sorted_subset.empty() || sorted_subset.back() != s.first) {
                    sorted_subset.push_back(s.first);
                }
                original_indices[s.second] = sorted_subset.size() - 1;
            }
        }
    }

    const std::vector<int>& extraction_subset() const {
        if (use_sorted_subset) {
            return sorted_subset;
        } else {
            return *original_subset;
        }
    }

    void fill_ranks(const double* ptr, RankedVector<double, int>& vec) const {
        if (use_sorted_subset) {
            size_t num = original_indices.size();
            for (size_t s = 0; s < num; ++s) {
                vec[s].first = ptr[original_indices[s]];
                vec[s].second = s;
            }
        } else {
            size_t num = original_subset->size();
            for (size_t s = 0; s < num; ++s, ++ptr) {
                vec[s].first = *ptr;
                vec[s].second = s;
            }
        }
        std::sort(vec.begin(), vec.end());
    }
};

template<typename Stat, typename Index>
void scaled_ranks(const RankedVector<Stat, Index>& collected, double* outgoing) { 
    // Computing tied ranks. 
    size_t cur_rank = 0;
    auto cIt = collected.begin();

    while (cIt != collected.end()) {
        auto copy = cIt;
        ++copy;
        double accumulated_rank = cur_rank;
        ++cur_rank;

        while (copy != collected.end() && copy->first == cIt->first) {
            accumulated_rank += cur_rank;
            ++cur_rank;
            ++copy;
        }

        double mean_rank= accumulated_rank / (copy - cIt);
        while (cIt!=copy) {
            outgoing[cIt->second] = mean_rank;
            ++cIt;
        }
    }

    // Mean-adjusting and converting to cosine values.
    double sum_squares = 0;
    size_t N = collected.size();
    const double center_rank = static_cast<double>(N - 1)/2; 
    for (size_t i = 0 ; i < N; ++i) {
        auto& o = outgoing[i];
        o -= center_rank;
        sum_squares += o*o;
    }

    // Special behaviour for no-variance cells; these are left as all-zero scaled ranks.
    sum_squares = std::max(sum_squares, 0.00000001);
    sum_squares = std::sqrt(sum_squares)*2;
    for (size_t i = 0; i < N; ++i) {
        outgoing[i] /= sum_squares;
    }

    return;
}

template<typename Stat, typename Index>
void subset_ranks(const RankedVector<Stat, Index>& x, RankedVector<Stat, Index>& output, const std::unordered_map<int, int>& subset) {
    for (size_t i = 0; i < x.size(); ++i) {
        auto it = subset.find(x[i].second);
        if (it != subset.end()) {
            output.emplace_back(x[i].first, it->second);
        }
    }
    return;
}

template<typename Stat, typename Index, typename Simple>
void simplify_ranks(const RankedVector<Stat, Index>& x, RankedVector<Simple, Index>& output) {
    if (x.size()) {
        Index counter = 0;
        auto last = x[0].first;
        for (const auto& r : x) {
            if (r.first != last) {
                ++counter;
                last = r.first;
            }
            output.emplace_back(counter, r.second);
        }
    }
    return;
}

}

#endif
