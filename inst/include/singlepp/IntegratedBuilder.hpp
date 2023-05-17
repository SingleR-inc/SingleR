#ifndef SINGLEPP_INTEGRATED_BUILDER_HPP
#define SINGLEPP_INTEGRATED_BUILDER_HPP

#include "macros.hpp"

#include "scaled_ranks.hpp"
#include "BasicBuilder.hpp"

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

/**
 * @file IntegratedBuilder.hpp
 *
 * @brief Prepare for integrated classification across references.
 */

namespace singlepp {

/**
 * @brief Reference datasets prepared for integrated classification.
 */
struct IntegratedReferences {
    /**
     * @return Number of reference datasets.
     * Each object corresponds to the reference used in an `IntegratedBuilder::add()` call, in the same order.
     */
    size_t num_references() const {
        return markers.size();
    }

    /**
     * @param r Reference dataset of interest.
     * @return Number of labels in this reference.
     */
    size_t num_labels(size_t r) const {
        return markers[r].size();
    }

    /**
     * @param r Reference dataset of interest.
     * @return Number of profiles in this reference.
     */
    size_t num_profiles(size_t r) const {
        size_t n = 0;
        for (const auto& ref : ranked[r]) {
            n += ref.size();
        }
        return n;
    }

    /**
     * @cond
     */
    std::vector<int> universe; // To be used by IntegratedScorer for indexed extraction.

    std::vector<bool> check_availability;
    std::vector<std::unordered_set<int> > available; // indices to 'universe'
    std::vector<std::vector<std::vector<int> > > markers; // indices to 'universe'
    std::vector<std::vector<std::vector<RankedVector<int, int> > > > ranked; // .second contains indices to 'universe'

    void resize(size_t n) {
        check_availability.resize(n);
        available.resize(n);
        markers.resize(n);
        ranked.resize(n);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Factory to prepare multiple references for integrated classification.
 *
 * For each reference dataset, we expect a `BasicBuilder::Prebuilt` or `BasicBuilder::PrebuiltIntersection` object,
 * as well as the original data structures (matrix, labels, etc.) used to construct that object.
 * These values are passed into `add()` to register that dataset, which can be repeated multiple times for different references.
 * Finally, calling `finish()` will return a vector of integrated references that can be used in `IntegratedScorer::run()`.
 * 
 * The preparation process mostly involves checking that the gene indices are consistent across references.
 * This is especially true when each reference contains a different set of features that must be intersected with the features in the test dataset.
 * See the documentation for `IntegratedScorer` for more details on the classification based on the integrated references.
 */
class IntegratedBuilder {
private:
    std::vector<const tatami::Matrix<double, int>*> stored_matrices;
    std::vector<const int*> stored_labels;
    IntegratedReferences references;
    std::vector<std::unordered_map<int, int> > gene_mapping;
    int nthreads = Defaults::num_threads;

public:
    /**
     * @brief Default parameters.
     */
    struct Defaults {
        /**
         * See `set_num_threads()` for details.
         */
        static constexpr int num_threads = 1;
    };

    /**
     * @param n Number of threads to use.
     *
     * @return A reference to this `IntegratedBuilder` object.
     */
    IntegratedBuilder& set_num_threads(int n = Defaults::num_threads) {
        nthreads = n;
        return *this;
    }

private:
    void add_internal_base(const tatami::Matrix<double, int>* ref, const int* labels) {
        stored_matrices.push_back(ref);
        stored_labels.push_back(labels);
        references.resize(stored_matrices.size());
        gene_mapping.resize(stored_matrices.size());
    }

    // Adding a reference without requiring any pruning of markers.
    template<class Subset>
    void add_internal_direct(const tatami::Matrix<double, int>* ref, const int* labels, const Markers& old_markers, const Subset& mat_subset) {
        add_internal_base(ref, labels);

        // Adding the markers for each label, indexed according to their
        // position in the test matrix. This assumes that 'mat_subset' is
        // appropriately specified to contain the test's row indices. 
        auto& new_markers = references.markers.back();
        new_markers.reserve(old_markers.size());

        for (size_t i = 0; i < old_markers.size(); ++i) {
            const auto& cur_old_markers = old_markers[i];

            std::unordered_set<int> unified;
            for (const auto& x : cur_old_markers) {
                unified.insert(x.begin(), x.end());
            }

            new_markers.emplace_back(unified.begin(), unified.end());

            if constexpr(!std::is_same<Subset, bool>::value) {
                auto& cur_new_markers = new_markers.back();
                for (auto& y : cur_new_markers) {
                    y = mat_subset[y];
                }
            }
        }
        return;
    }

    // Adding a reference with different features from the test, requiring some
    // pruning to remove markers that are not present in the intersection.
    template<class Subset>
    void add_internal_intersect(
        const std::vector<std::pair<int, int> >& intersection,
        const tatami::Matrix<double, int>* ref, 
        const int* labels, 
        const Markers& old_markers,
        const Subset& ref_subset)
    {
        add_internal_base(ref, labels);

        // Manually constructing the markers. This involves (i) pruning out the
        // markers that aren't present in the intersection, and (ii) updating 
        // their indices so that they point to rows of 'mat', not 'ref'.
        std::unordered_map<int, int> reverse_map;
        for (const auto& i : intersection) {
            reverse_map[i.second] = i.first;
        }

        auto subindex = [&](int i) -> int {
            if constexpr(!std::is_same<Subset, bool>::value) {
                return ref_subset[i];
            } else {
                return i;
            }
        };

        auto& new_markers = references.markers.back();
        new_markers.resize(old_markers.size());

        for (size_t i = 0; i < old_markers.size(); ++i) {
            const auto& cur_old_markers = old_markers[i];

            std::unordered_set<int> unified;
            for (const auto& x : cur_old_markers) {
                unified.insert(x.begin(), x.end());
            }

            auto& cur_new_markers = new_markers[i];
            cur_new_markers.reserve(unified.size());

            for (auto y : unified) {
                auto it = reverse_map.find(subindex(y));
                if (it != reverse_map.end()) {
                    cur_new_markers.push_back(it->second);
                }
            }
        }

        // Constructing the mapping of mat's rows to the reference rows.
        references.check_availability.back() = true;
        auto& mapping = gene_mapping.back();
        for (const auto& i : intersection) {
            mapping[i.first] = i.second;
        }
        return;
    }

public:
    /**
     * @param ref Matrix containing the reference expression values.
     * Rows are features and columns are reference samples.
     * The number and identity of features should be identical to the test dataset to be classified in `IntegratedScorer`.
     * @param[in] labels Pointer to an array of label assignments.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param markers A vector of vectors of ranked marker genes for each pairwise comparison between labels in `ref`, see `Markers` for more details.
     *
     * @return The reference dataset is registered for later use in `finish()`.
     *
     * `ref` and `labels` are expected to remain valid until `finish()` is called.
     */
    void add(const tatami::Matrix<double, int>* ref, const int* labels, const Markers& markers) {
        add_internal_direct(ref, labels, markers, false);    
    }

    /**
     * @tparam Id Type of the gene identifier for each row.
     *
     * @param mat_nrow Number of rows (genes) in the test dataset.
     * @param[in] mat_id Pointer to an array of identifiers of length equal to `mat_nrow`.
     * This should contain a unique identifier for each row of `mat` (typically a gene name or index).
     * If any duplicate IDs are present, only the first occurrence is used.
     * @param ref An expression matrix for the reference expression profiles, where rows are genes and columns are cells.
     * This should have non-zero columns.
     * @param[in] ref_id Pointer to an array of identifiers of length equal to the number of rows of any `ref`.
     * This should contain a unique identifier for each row in `ref`, and should be comparable to `mat_id`.
     * If any duplicate IDs are present, only the first occurrence is used.
     * @param[in] labels An array of length equal to the number of columns of `ref`, containing the label for each sample.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param markers A vector of vectors of ranked marker genes for each pairwise comparison between labels in `ref`, see `Markers` for more details.
     *
     * @return The reference dataset is registered for later use in `finish()`.
     *
     * `ref` and `labels` are expected to remain valid until `finish()` is called.
     * `mat_id` and `mat_nrow` should also be constant for all invocations to `add()`.
     */
    template<typename Id>
    void add(size_t mat_nrow,
        const Id* mat_id,
        const tatami::Matrix<double, int>* ref, 
        const Id* ref_id,
        const int* labels, 
        const Markers& markers)
    {
        auto intersection = intersect_features(mat_nrow, mat_id, ref->nrow(), ref_id);
        add_internal_intersect(intersection, ref, labels, markers, false);
    }

public:
    /**
     * @param ref Matrix containing the reference expression values.
     * Rows are features and columns are reference samples.
     * The number and identity of features should be identical to the test dataset to be classified in `IntegratedScorer`.
     * @param[in] labels Pointer to an array of label assignments.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param built The built reference created by running `BasicBuilder::run()` on `ref` and `labels`.
     *
     * @return The reference dataset is registered for later use in `finish()`.
     *
     * `ref` and `labels` are expected to remain valid until `finish()` is called.
     */
    void add(const tatami::Matrix<double, int>* ref, const int* labels, const BasicBuilder::Prebuilt& built) {
        add_internal_direct(ref, labels, built.markers, built.subset);
        return;
    }

    /**
     * @param intersection Vector defining the intersection of features betweent the test and reference datasets.
     * Each entry is a pair where the first element is the row index in the test matrix,
     * and the second element is the row index for the corresponding feature in the reference matrix.
     * Each row index for either matrix should occur no more than once in `intersection`.
     * @param ref An expression matrix for the reference expression profiles, where rows are genes and columns are cells.
     * This should have non-zero columns.
     * @param[in] labels An array of length equal to the number of columns of `ref`, containing the label for each sample.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param built The built reference created by running `BasicBuilder::run()` on all preceding arguments.
     *
     * @return The reference dataset is registered for later use in `finish()`.
     *
     * `ref` and `labels` are expected to remain valid until `finish()` is called.
     * `mat_id` and `mat_nrow` should also be constant for all invocations to `add()`.
     */
    void add(const std::vector<std::pair<int, int> >& intersection,
        const tatami::Matrix<double, int>* ref, 
        const int* labels, 
        const BasicBuilder::PrebuiltIntersection& built) 
    {
        add_internal_direct(ref, labels, built.markers, built.mat_subset);
        references.check_availability.back() = true;

        // Constructing the mapping of mat's rows to the reference rows.
        auto& mapping = gene_mapping.back();
        for (const auto& i : intersection) {
            mapping[i.first] = i.second;
        }
        return;
    }

    /**
     * @tparam Id Type of the gene identifier for each row.
     *
     * @param mat_nrow Number of rows (genes) in the test dataset.
     * @param[in] mat_id Pointer to an array of identifiers of length equal to `mat_nrow`.
     * This should contain a unique identifier for each row of `mat` (typically a gene name or index).
     * If any duplicate IDs are present, only the first occurrence is used.
     * @param ref An expression matrix for the reference expression profiles, where rows are genes and columns are cells.
     * This should have non-zero columns.
     * @param[in] ref_id Pointer to an array of identifiers of length equal to the number of rows of any `ref`.
     * This should contain a unique identifier for each row in `ref`, and should be comparable to `mat_id`.
     * If any duplicate IDs are present, only the first occurrence is used.
     * @param[in] labels An array of length equal to the number of columns of `ref`, containing the label for each sample.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param built The built reference created by running `BasicBuilder::run()` on all preceding arguments.
     *
     * @return The reference dataset is registered for later use in `finish()`.
     *
     * `ref` and `labels` are expected to remain valid until `finish()` is called.
     * `mat_id` and `mat_nrow` should also be constant for all invocations to `add()`.
     */
    template<typename Id>
    void add(size_t mat_nrow,
        const Id* mat_id,
        const tatami::Matrix<double, int>* ref, 
        const Id* ref_id,
        const int* labels, 
        const BasicBuilder::PrebuiltIntersection& built) 
    {
        auto intersection = intersect_features(mat_nrow, mat_id, ref->nrow(), ref_id);
        add(intersection, ref, labels, built);
        return;
    }

    /**
     * @param intersection Vector defining the intersection of features betweent the test and reference datasets.
     * Each entry is a pair where the first element is the row index in the test matrix,
     * and the second element is the row index for the corresponding feature in the reference matrix.
     * Each row index for either matrix should occur no more than once in `intersection`.
     * @param ref An expression matrix for the reference expression profiles, where rows are genes and columns are cells.
     * This should have non-zero columns.
     * @param[in] labels An array of length equal to the number of columns of `ref`, containing the label for each sample.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param built The built reference created by running `BasicBuilder::run()` on `ref` and `labels`.
     *
     * @return The reference dataset is registered for later use in `finish()`.
     *
     * `ref` and `labels` are expected to remain valid until `finish()` is called.
     * `mat_id` and `mat_nrow` should also be constant for all invocations to `add()`.
     */
    void add(const std::vector<std::pair<int, int> >& intersection,
        const tatami::Matrix<double, int>* ref, 
        const int* labels, 
        const BasicBuilder::Prebuilt& built) 
    {
        add_internal_intersect(intersection, ref, labels, built.markers, built.subset);
    }

    /**
     * @tparam Id Type of the gene identifier for each row.
     *
     * @param mat_nrow Number of rows (genes) in the test dataset.
     * @param[in] mat_id Pointer to an array of identifiers of length equal to `mat_nrow`.
     * This should contain a unique identifier for each row of `mat` (typically a gene name or index).
     * If any duplicate IDs are present, only the first occurrence is used.
     * @param ref An expression matrix for the reference expression profiles, where rows are genes and columns are cells.
     * This should have non-zero columns.
     * @param[in] ref_id Pointer to an array of identifiers of length equal to the number of rows of any `ref`.
     * This should contain a unique identifier for each row in `ref`, and should be comparable to `mat_id`.
     * If any duplicate IDs are present, only the first occurrence is used.
     * @param[in] labels An array of length equal to the number of columns of `ref`, containing the label for each sample.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param built The built reference created by running `BasicBuilder::run()` on `ref` and `labels`.
     *
     * @return The reference dataset is registered for later use in `finish()`.
     *
     * `ref` and `labels` are expected to remain valid until `finish()` is called.
     * `mat_id` and `mat_nrow` should also be constant for all invocations to `add()`.
     */
    template<typename Id>
    void add(size_t mat_nrow,
        const Id* mat_id,
        const tatami::Matrix<double, int>* ref, 
        const Id* ref_id,
        const int* labels, 
        const BasicBuilder::Prebuilt& built) 
    {
        auto intersection = intersect_features(mat_nrow, mat_id, ref->nrow(), ref_id);
        add(intersection, ref, labels, built);
    }

private:
    /* Here, we've split out some of the functions for easier reading.
     * Otherwise everything would land in a single mega-function.
     */

    static void fill_ranks(
        const tatami::Matrix<double, int>* curmat, 
        const int* curlab, 
        const std::vector<int>& in_use, 
        const std::vector<int>& positions, 
        std::vector<std::vector<RankedVector<int, int> > >& cur_ranked,
        int nthreads) 
    {
        // If we don't need to check availability, this implies that 
        // the reference has 1:1 feature mapping to the test dataset.
        // In that case, we can proceed quite simply.
        tatami::parallelize([&](int, int start, int len) -> void {
            RankedVector<double, int> tmp_ranked;
            tmp_ranked.reserve(in_use.size());

            // 'in_use' is guaranteed to be sorted and unique, see its derivation in finish().
            // This means we can directly use it for indexed extraction.
            auto wrk = tatami::consecutive_extractor<false, false>(curmat, start, len, in_use); 
            std::vector<double> buffer(wrk->index_length);

            for (int c = start, end = start + len; c < end; ++c) {
                auto ptr = wrk->fetch(c, buffer.data());

                tmp_ranked.clear();
                for (int i = 0, end = in_use.size(); i < end; ++i, ++ptr) {
                    tmp_ranked.emplace_back(*ptr, i);
                }
                std::sort(tmp_ranked.begin(), tmp_ranked.end());

                auto& final_ranked = cur_ranked[curlab[c]][positions[c]];
                simplify_ranks(tmp_ranked, final_ranked);
            }
        }, curmat->ncol(), nthreads);
    }

    static void fill_ranks(
        const tatami::Matrix<double, int>* curmat, 
        const int* curlab, 
        const std::vector<int>& in_use, 
        const std::vector<int>& positions,
        const std::unordered_map<int, int>& cur_mapping,
        std::unordered_set<int>& cur_available,
        std::vector<std::vector<RankedVector<int, int> > >& cur_ranked,
        int nthreads)
    {
        // If we need to check availability, then we need to check
        // the mapping of test genes to row indices of the reference.
        std::vector<std::pair<int, int> > remapping; 
        remapping.reserve(in_use.size());

        for (int i = 0, end = in_use.size(); i < end; ++i) {
            auto it = cur_mapping.find(in_use[i]);
            if (it != cur_mapping.end()) {
                remapping.emplace_back(it->second, i); // using 'i' instead of 'in_use[i]', as we want to work with indices into 'in_use', not the values of 'in_use' themselves.
                cur_available.insert(i);
            }
        }

        std::sort(remapping.begin(), remapping.end());

        // This section is just to enable indexed extraction by tatami.
        // There's no need to consider duplicates among the
        // 'remapping[i].first', 'cur_mapping->second' is guaranteed to be
        // unique as a consequence of how intersect_features works.
        std::vector<int> remapped_in_use; 
        remapped_in_use.reserve(remapping.size());
        for (const auto& p : remapping) {
            remapped_in_use.push_back(p.first);
        }

        tatami::parallelize([&](int, int start, int len) -> void {
            RankedVector<double, int> tmp_ranked;
            tmp_ranked.reserve(remapped_in_use.size());
            std::vector<double> buffer(remapped_in_use.size());
            auto wrk = tatami::consecutive_extractor<false, false>(curmat, start, len, remapped_in_use);

            for (size_t c = start, end = start + len; c < end; ++c) {
                auto ptr = wrk->fetch(c, buffer.data());

                tmp_ranked.clear();
                for (const auto& p : remapping) {
                    tmp_ranked.emplace_back(*ptr, p.second); // remember, 'p.second' corresponds to indices into 'in_use'.
                    ++ptr;
                }
                std::sort(tmp_ranked.begin(), tmp_ranked.end());

                auto& final_ranked = cur_ranked[curlab[c]][positions[c]];
                simplify_ranks(tmp_ranked, final_ranked);
            }
        }, curmat->ncol(), nthreads);
    }

public:
    /**
     * @return The set of integrated references, for use in `IntegratedScorer`.
     *
     * This function should only be called once, after all reference datasets have been registered with `add()`.
     * Any further invocations of this function will not be valid.
     */
    IntegratedReferences finish() {
        // Identify the union of all marker genes.
        auto& in_use = references.universe;
        std::unordered_map<int, int> remap_to_universe;
        {
            std::unordered_set<int> in_use_tmp;
            for (const auto& refmarkers : references.markers) {
                for (const auto& mrk : refmarkers) {
                    in_use_tmp.insert(mrk.begin(), mrk.end());
                }
            }

            in_use.insert(in_use.end(), in_use_tmp.begin(), in_use_tmp.end());
            std::sort(in_use.begin(), in_use.end());

            for (int i = 0, end = in_use.size(); i < end; ++i) {
                remap_to_universe[in_use[i]] = i;
            }
        }

        for (size_t r = 0, end = references.num_references(); r < end; ++r) {
            auto curlab = stored_labels[r];
            auto curmat = stored_matrices[r];

            // Reindexing the markers to point to the universe.
            auto& curmarkers = references.markers[r];
            for (auto& outer : curmarkers) {
                for (auto& x : outer) {
                    x = remap_to_universe.find(x)->second;
                }
            }

            auto& cur_ranked = references.ranked[r];
            std::vector<int> positions;
            {
                size_t nlabels = curmarkers.size();
                size_t NC = curmat->ncol();
                positions.reserve(NC);                

                std::vector<int> samples_per_label(nlabels);
                for (size_t c = 0; c < NC; ++c) {
                    auto& pos = samples_per_label[curlab[c]];
                    positions.push_back(pos);
                    ++pos;
                }

                cur_ranked.resize(nlabels);
                for (size_t l = 0; l < nlabels; ++l) {
                    cur_ranked[l].resize(samples_per_label[l]);
                }
            }

            // Finally filling the rankings.
            if (!references.check_availability[r]) {
                fill_ranks(curmat, curlab, in_use, positions, references.ranked[r], nthreads);
            } else {
                fill_ranks(curmat, curlab, in_use, positions, gene_mapping[r], references.available[r], references.ranked[r], nthreads);
            }
        }

        return std::move(references);
    }
};

}

#endif
