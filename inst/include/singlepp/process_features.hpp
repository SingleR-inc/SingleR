#ifndef SINGLEPP_PROCESS_FEATURES_HPP
#define SINGLEPP_PROCESS_FEATURES_HPP

#include "macros.hpp"

#include "Markers.hpp"

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

namespace singlepp {

typedef std::vector<std::pair<int, int> > Intersection;

template<typename Id>
Intersection intersect_features(size_t mat_n, const Id* mat_id, size_t ref_n, const Id* ref_id) {
    std::unordered_map<Id, std::pair<int, int> > intersection;
    intersection.reserve(mat_n);
    for (size_t i = 0; i < mat_n; ++i) {
        if (intersection.find(mat_id[i]) == intersection.end()) { // only using the first occurrence of each ID in mat_id.
            intersection[mat_id[i]] = std::make_pair<int, int>(i, -1);
        }
    }

    for (size_t i = 0; i < ref_n; ++i) {
        auto it = intersection.find(ref_id[i]);
        if (it != intersection.end()) { // only using the first occurrence of each ID in ref_id.
            auto& target = (it->second).second;
            if (target < 0) {
                target = i;
            }
        }
    }

    Intersection pairings;
    pairings.reserve(intersection.size());
    for (const auto& x : intersection) {
        if (x.second.second >= 0) {
            pairings.push_back(x.second);
        }
    }

    return pairings;
}

inline void subset_markers(Intersection& intersection, Markers& markers, int top) {
    std::unordered_set<int> available;
    available.reserve(intersection.size());
    for (const auto& in : intersection) {
        available.insert(in.second);
    }

    // Figuring out the top markers to retain, that are _also_ in the intersection.
    std::unordered_set<int> all_markers;
    all_markers.reserve(intersection.size());

    for (size_t i = 0; i < markers.size(); ++i) {
        for (size_t j = 0; j < markers[i].size(); ++j) {
            auto& current = markers[i][j];

            std::vector<int> replacement;
            size_t upper_bound = static_cast<size_t>(top >= 0 ? top : -1); // in effect, no upper bound if top = -1.
            replacement.reserve(top >= 0 ? static_cast<size_t>(top) : current.size());

            for (size_t k = 0; k < current.size() && replacement.size() < upper_bound; ++k) {
                if (available.find(current[k]) != available.end()) {
                    all_markers.insert(current[k]);
                    replacement.push_back(current[k]);
                }
            }
            current.swap(replacement);
        }
    }

    // Subsetting the intersection down to the chosen set of markers.
    std::unordered_map<int, int> mapping;
    mapping.reserve(intersection.size());
    {
        size_t counter = 0;
        for (size_t i = 0; i < intersection.size(); ++i) {
            if (all_markers.find(intersection[i].second) != all_markers.end()) {
                intersection[counter] = intersection[i];
                mapping[intersection[i].second] = counter;
                ++counter;
            }
        }
        intersection.resize(counter);
    }

    // Reindexing the markers.
    for (size_t i = 0; i < markers.size(); ++i) {
        for (size_t j = 0; j < markers[i].size(); ++j) {
            auto& current = markers[i][j];
            for (size_t k = 0; k < current.size(); ++k) {
                auto it = mapping.find(current[k]);
                current[k] = it->second;
            }
        }
    }

    return;
}

// Use this method when the feature spaces are already identical.
inline std::vector<int> subset_markers(Markers& markers, int top) {
    std::unordered_set<int> available;
    for (size_t i = 0; i < markers.size(); ++i) {
        for (size_t j = 0; j < markers[i].size(); ++j) {
            auto& current = markers[i][j];
            if (top >= 0) {
                current.resize(std::min(current.size(), static_cast<size_t>(top)));
            }
            available.insert(current.begin(), current.end());
        }
    }

    std::vector<int> subset(available.begin(), available.end());
    std::sort(subset.begin(), subset.end());

    std::unordered_map<int, int> mapping;
    mapping.reserve(subset.size());
    for (size_t i = 0; i < subset.size(); ++i) {
        mapping[subset[i]] = i;
    }

    // Reindexing the markers.
    for (size_t i = 0; i < markers.size(); ++i) {
        for (size_t j = 0; j < markers[i].size(); ++j) {
            auto& current = markers[i][j];
            for (size_t k = 0; k < current.size(); ++k) {
                auto it = mapping.find(current[k]);
                current[k] = it->second;
            }
        }
    }

    return subset;
}

inline std::pair<std::vector<int>, std::vector<int> > unzip(const Intersection& intersection) {
    std::vector<int> left(intersection.size()), right(intersection.size());
    for (size_t i = 0; i < intersection.size(); ++i) {
        left[i] = intersection[i].first;
        right[i] = intersection[i].second;
    }
    return std::make_pair(std::move(left), std::move(right));
}

}

#endif
