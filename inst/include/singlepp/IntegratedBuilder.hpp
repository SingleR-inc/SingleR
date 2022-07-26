#ifndef SINGLEPP_INTEGRATED_BUILDER_HPP
#define SINGLEPP_INTEGRATED_BUILDER_HPP

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
 * @brief Single reference dataset prepared for integrated classification.
 */
struct IntegratedReference {
    /**
     * @return Number of labels in this reference.
     */
    size_t num_labels() const {
        return markers.size();
    }

    /**
     * @return Number of profiles in this reference.
     */
    size_t num_profiles() const {
        size_t n = 0;
        for (const auto& ref : ranked) {
            n += ref.size();
        }
        return n;
    }

    /**
     * @cond
     */
    bool check_availability = false;
    std::unordered_set<int> available;
    std::vector<std::vector<int> > markers;
    std::vector<std::vector<RankedVector<int, int> > > ranked;
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
    std::vector<IntegratedReference> references;
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
    void add_internal(const tatami::Matrix<double, int>* ref, const int* labels) {
        stored_matrices.push_back(ref);
        stored_labels.push_back(labels);
        references.resize(references.size() + 1);
        gene_mapping.resize(gene_mapping.size() + 1); 
    }

    template<class Subset>
    void add_internal(const tatami::Matrix<double, int>* ref, const int* labels, const Markers& old_markers, const Subset& subset) {
        add_internal(ref, labels);

        // Adding the markers for each label, indexed according to their
        // position in the test matrix. This assumes that 'subset' is
        // appropriately specified to contain the test's row indices. Note that
        // these are positions in the test, not the ref; IntegratedScorer will
        // use the gene_mapping to convert them to row indices of 'ref'.
        auto& new_markers = references.back().markers;
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
                    y = subset[y];
                }
            }
        }
        return;
    }

    template<typename Id, class Subset>
    void add_internal(size_t mat_nrow,
        const Id* mat_id,
        const tatami::Matrix<double, int>* ref, 
        const Id* ref_id,
        const int* labels, 
        const Markers& old_markers,
        const Subset& subset)
    {
        add_internal(ref, labels);

        auto intersection = intersect_features(mat_nrow, mat_id, ref->nrow(), ref_id);

        // Manually constructing the markers. This involves (i) pruning out the
        // markers that aren't present in the intersection, and (ii) updating 
        // their indices so that they point to rows of 'mat', not 'ref'.
        std::unordered_map<int, int> reverse_map;
        for (const auto& i : intersection) {
            reverse_map[i.second] = i.first;
        }

        auto subindex = [&](int i) -> int {
            if constexpr(!std::is_same<Subset, bool>::value) {
                return subset[i];
            } else {
                return i;
            }
        };

        auto& new_markers = references.back().markers;
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
        references.back().check_availability = true;
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
        add_internal(ref, labels, markers, false);    
    }

    /**
     * @tparam Id Type of the gene identifier for each row.
     *
     * @param mat_nrow Number of rows (genes) in the test dataset.
     * @param[in] mat_id Pointer to an array of identifiers of length equal to `mat_nrow`.
     * This should contain a unique identifier for each row of `mat` (typically a gene name or index).
     * @param ref An expression matrix for the reference expression profiles, where rows are genes and columns are cells.
     * This should have non-zero columns.
     * @param[in] ref_id Pointer to an array of identifiers of length equal to the number of rows of any `ref`.
     * This should contain a unique identifier for each row in `ref`, and should be comparable to `mat_id`.
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
        add_internal(mat_nrow, mat_id, ref, ref_id, labels, markers, false);
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
        add_internal(ref, labels, built.markers, built.subset);
        return;
    }

    /**
     * @tparam Id Type of the gene identifier for each row.
     *
     * @param mat_nrow Number of rows (genes) in the test dataset.
     * @param[in] mat_id Pointer to an array of identifiers of length equal to `mat_nrow`.
     * This should contain a unique identifier for each row of `mat` (typically a gene name or index).
     * @param ref An expression matrix for the reference expression profiles, where rows are genes and columns are cells.
     * This should have non-zero columns.
     * @param[in] ref_id Pointer to an array of identifiers of length equal to the number of rows of any `ref`.
     * This should contain a unique identifier for each row in `ref`, and should be comparable to `mat_id`.
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
        add_internal(ref, labels, built.markers, built.mat_subset);
        references.back().check_availability = true;

        // Constructing the mapping of mat's rows to the reference rows.
        auto intersection = intersect_features(mat_nrow, mat_id, ref->nrow(), ref_id);
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
     * @param ref An expression matrix for the reference expression profiles, where rows are genes and columns are cells.
     * This should have non-zero columns.
     * @param[in] ref_id Pointer to an array of identifiers of length equal to the number of rows of any `ref`.
     * This should contain a unique identifier for each row in `ref`, and should be comparable to `mat_id`.
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
        add_internal(mat_nrow, mat_id, ref, ref_id, labels, built.markers, built.subset);
    }

public:
    /**
     * @return A vector of `IntegratedReference` objects.
     * Each object corresponds to the reference used in an `add()` call, in the same order.
     *
     * This function should only be called once, after all reference datasets have been registered with `add()`.
     * Any further invocations of this function will not be valid.
     */
    std::vector<IntegratedReference> finish() {
        /**
         * @cond
         */
        // Identify the global set of all genes that will be in use here.
        std::unordered_set<int> in_use_tmp;
        for (const auto& ref : references) {
            for (const auto& mrk : ref.markers) {
                in_use_tmp.insert(mrk.begin(), mrk.end());
            }
        }

        std::vector<int> in_use(in_use_tmp.begin(), in_use_tmp.end());
        std::sort(in_use.begin(), in_use.end());

        for (size_t r = 0; r < references.size(); ++r) {
            auto& curref = references[r];
            auto curlab = stored_labels[r];
            auto curmat = stored_matrices[r];

            size_t NR = curmat->nrow();
            size_t NC = curmat->ncol();
            size_t nlabels = curref.markers.size();

            std::vector<int> positions(NC);
            std::vector<int> samples_per_label(nlabels);
            for (size_t c = 0; c < NC; ++c) {
                auto& pos = samples_per_label[curlab[c]];
                positions[c] = pos;
                ++pos;
            }

            curref.ranked.resize(nlabels);
            for (size_t l = 0; l < nlabels; ++l) {
                curref.ranked[l].resize(samples_per_label[l]);
            }

            if (!curref.check_availability) {
                // If we don't need to check availability, this implies that 
                // the reference has 1:1 feature mapping to the test dataset.
                // In that case, we can proceed quite simply.
                size_t first = 0, last = 0;
                if (in_use.size()) {
                    first = in_use.front();
                    last = in_use.back() + 1;
                }

#ifndef SINGLEPP_CUSTOM_PARALLEL
                #pragma omp parallel num_threads(nthreads)
                {
#else
                SINGLEPP_CUSTOM_PARALLEL(NC, [&](size_t start, size_t end) -> void {
#endif

                    RankedVector<double, int> tmp_ranked;
                    tmp_ranked.reserve(in_use.size());
                    std::vector<double> buffer(NR);
                    auto wrk = curmat->new_workspace(false);

#ifndef SINGLEPP_CUSTOM_PARALLEL
                    #pragma omp for
                    for (size_t c = 0; c < NC; ++c) {
#else
                    for (size_t c = start; c < end; ++c) {
#endif

                        auto ptr = curmat->column(c, buffer.data(), first, last, wrk.get());

                        tmp_ranked.clear();
                        for (auto u : in_use) {
                            tmp_ranked.emplace_back(ptr[u - first], u);
                        }
                        std::sort(tmp_ranked.begin(), tmp_ranked.end());

                        auto& final_ranked = curref.ranked[curlab[c]][positions[c]];
                        simplify_ranks(tmp_ranked, final_ranked);
                    }

#ifndef SINGLEPP_CUSTOM_PARALLEL
                }
#else
                }, nthreads);
#endif

            } else {
                // If we do need to check availability, then we need to check
                // the mapping of test genes to their reference row indices.
                const auto& cur_mapping = gene_mapping[r];
                auto& cur_available = curref.available;
                std::unordered_map<int, int> remapping;
                remapping.reserve(in_use.size());
                size_t first = NR, last = 0;

                for (auto u : in_use) {
                    auto it = cur_mapping.find(u);
                    if (it == cur_mapping.end()) {
                        continue;
                    }

                    remapping[u] = it->second;
                    cur_available.insert(u);
                    if (it->second < first) {
                        first = it->second;
                    }
                    if (it->second > last) {
                        last = it->second;
                    }
                }
                last = std::max(last + 1, first);

#ifndef SINGLEPP_CUSTOM_PARALLEL
                #pragma omp parallel num_threads(nthreads)
                {
#else
                SINGLEPP_CUSTOM_PARALLEL(NC, [&](size_t start, size_t end) -> void {
#endif

                    RankedVector<double, int> tmp_ranked;
                    tmp_ranked.reserve(in_use.size());
                    std::vector<double> buffer(NR);
                    auto wrk = curmat->new_workspace(false);

#ifndef SINGLEPP_CUSTOM_PARALLEL
                    #pragma omp for
                    for (size_t c = 0; c < NC; ++c) {
#else
                    for (size_t c = start; c < end; ++c) {
#endif

                        auto ptr = curmat->column(c, buffer.data(), first, last, wrk.get());

                        tmp_ranked.clear();
                        for (auto p : remapping) {
                            tmp_ranked.emplace_back(ptr[p.second - first], p.first);
                        }
                        std::sort(tmp_ranked.begin(), tmp_ranked.end());

                        auto& final_ranked = curref.ranked[curlab[c]][positions[c]];
                        simplify_ranks(tmp_ranked, final_ranked);
                    }

#ifndef SINGLEPP_CUSTOM_PARALLEL
                }
#else
                }, nthreads);
#endif
            }
        }
        /**
         * @endcond
         */

        return std::move(references);
    }
};

}

#endif
