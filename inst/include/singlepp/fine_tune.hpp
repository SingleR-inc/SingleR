#ifndef TATAMI_FINE_TUNE_HPP
#define TATAMI_FINE_TUNE_HPP

#include "macros.hpp"

#include "scaled_ranks.hpp"
#include "process_features.hpp"
#include "compute_scores.hpp"
#include "build_indices.hpp"

#include <vector>
#include <algorithm>
#include <unordered_map>

namespace singlepp {

inline std::pair<int, double> fill_labels_in_use(const std::vector<double>& scores, double threshold, std::vector<int>& in_use) {
    auto it = std::max_element(scores.begin(), scores.end());
    int best_label = it - scores.begin();
    double max_score = *it;

    in_use.clear();
    constexpr double DUMMY = -1000;
    double next_score = DUMMY;
    const double bound = max_score - threshold;

    for (size_t i = 0; i < scores.size(); ++i) {
        const auto& val = scores[i];
        if (val >= bound) {
            in_use.push_back(i);
        }
        if (i != best_label && next_score < val) {
            next_score = val;
        }
    }

    return std::make_pair(best_label, max_score - next_score); 
}

inline std::pair<int, double> replace_labels_in_use(const std::vector<double>& scores, double threshold, std::vector<int>& in_use) {
    auto it = std::max_element(scores.begin(), scores.end());
    int best_index = it - scores.begin();
    double max_score = *it;

    int best_label = in_use[best_index];
    size_t counter = 0;

    constexpr double DUMMY = -1000;
    double next_score = DUMMY;
    const double bound = max_score - threshold;

    for (size_t i = 0; i < scores.size(); ++i) {
        const auto& val = scores[i];
        if (val >= bound) {
            in_use[counter] = in_use[i];
            ++counter;
        }
        if (i != best_index && next_score < val) {
            next_score = val;
        }
    }

    in_use.resize(counter);
    return std::make_pair(best_label, max_score - next_score); 
}

class FineTuner {
    std::vector<int> labels_in_use;

    std::unordered_map<int, int> gene_subset;

    std::vector<double> scaled_left, scaled_right;

    std::vector<double> all_correlations;

    RankedVector<double, int> input_sub;

    RankedVector<int, int> ref_sub;

public:
    template<bool test = false>
    std::pair<int, double> run(
        const RankedVector<double, int>& input, 
        const std::vector<Reference>& ref,
        const Markers& markers,
        std::vector<double>& scores,
        double quantile,
        double threshold)
    {
        if (scores.size() <= 1) {
            return std::make_pair(0, std::numeric_limits<double>::quiet_NaN());
        } 

        auto candidate = fill_labels_in_use(scores, threshold, labels_in_use);

        // If there's only one top label, we don't need to do anything else.
        if (labels_in_use.size() == 1) {
            return candidate;
        } 

        // We also give up if every label is in range, because any subsequent
        // calculations would use all markers and just give the same result.
        // The 'test' parameter allows us to skip this bypass for testing.
        if constexpr(!test) {
            if (labels_in_use.size() == ref.size()) {
                return candidate;
            }
        }

        while (labels_in_use.size() > 1) {
            gene_subset.clear();
            gene_subset.reserve(input.size()); // known maximum.
            for (auto l : labels_in_use) {
                for (auto l2 : labels_in_use){ 
                    const auto& current = markers[l][l2];
                    for (auto c : current) {
                        if (gene_subset.find(c) == gene_subset.end()) {
                            const int counter = gene_subset.size();
                            gene_subset[c] = counter;
                        }
                    }
                }
            }

            input_sub.clear();
            subset_ranks(input, input_sub, gene_subset);
            scaled_left.resize(gene_subset.size());
            scaled_ranks(input_sub, scaled_left.data());

            scaled_right.resize(gene_subset.size());
            scores.clear();

            for (size_t i = 0; i < labels_in_use.size(); ++i) {
                auto curlab = labels_in_use[i];

                all_correlations.clear();
                const auto& curref = ref[curlab];
                size_t NR = curref.index->ndim();
                size_t NC = curref.index->nobs();

                for (size_t c = 0; c < NC; ++c) {
                    // Technically we could be faster if we remembered the
                    // subset from the previous fine-tuning iteration, but this
                    // requires us to (possibly) make a copy of the entire
                    // reference set; we can't afford to do this in each thread.
                    ref_sub.clear();
                    subset_ranks(curref.ranked[c], ref_sub, gene_subset);
                    scaled_ranks(ref_sub, scaled_right.data());

                    double cor = distance_to_correlation(scaled_left.size(), scaled_left, scaled_right);
                    all_correlations.push_back(cor);
                }

                double score = correlations_to_scores(all_correlations, quantile);
                scores.push_back(score);
            }

            candidate = replace_labels_in_use(scores, threshold, labels_in_use); 
            if (labels_in_use.size() == scores.size()) { // i.e., unchanged.
                break;
            }
        }

        return candidate;
    }
};

}

#endif
