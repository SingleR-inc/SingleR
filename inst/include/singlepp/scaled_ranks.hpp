#ifndef SINGLEPP_SCALED_RANKS_HPP
#define SINGLEPP_SCALED_RANKS_HPP

#include <algorithm>
#include <vector>
#include <cmath>
#include <unordered_map>

namespace singlepp {

template<typename Stat, typename Index>
using RankedVector = std::vector<std::pair<Stat, Index> >;

template<typename Stat, typename Index>
void fill_ranks(size_t num_subset, const int* subset, const Stat* ptr, RankedVector<Stat, Index>& vec, int offset = 0) {
    for (size_t s = 0; s < num_subset; ++s, ++subset) {
        vec[s].first = ptr[*subset - offset];
        vec[s].second = s;
    }
    std::sort(vec.begin(), vec.end());
    return;
}

template<typename Stat, typename Index>
void fill_ranks(const std::vector<int>& subset, const Stat* ptr, RankedVector<Stat, Index>& vec, int offset = 0) {
    fill_ranks(subset.size(), subset.data(), ptr, vec, offset);
    return;
}

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
