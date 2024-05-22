#ifndef SINGLEPP_BASIC_BUILDER_HPP
#define SINGLEPP_BASIC_BUILDER_HPP

#include "macros.hpp"

#include "knncolle/knncolle.hpp"
#include "tatami/tatami.hpp"

#include "build_indices.hpp"
#include "process_features.hpp"

#include <vector>

/**
 * @file BasicBuilder.hpp
 *
 * @brief Defines the `BasicBuilder` class.
 */

namespace singlepp {

/**
 * @brief Construct a single pre-built reference for further classification.
 *
 * This class prepares a single labelled reference dataset for classification of a test dataset.
 * We select the top marker genes to use for each pairwise comparison between labels,
 * and we use them to pre-compute the ranks in advance of computing the Spearman correlations.
 * The pre-built reference can then be used in `BasicScorer::run()`. 
 */
class BasicBuilder {
public:
    /**
     * @brief Default parameters for reference building.
     */
    struct Defaults {
        /**
         * See `set_top()` for details.
         */
        static constexpr int top = -1;

        /**
         * See `set_approximate()` for details.
         */
        static constexpr bool approximate = false;

        /**
         * See `set_num_threads()` for details.
         */
        static constexpr int num_threads = 1;
    };

private:
    int top = Defaults::top;
    bool approximate = Defaults::approximate;
    int nthreads = Defaults::num_threads;

public:
    /**
     * @param t Number of top markers to use from each pairwise comparison between labels, see `Classifier::set_top()`.
     *
     * @return A reference to this `BasicBuilder` object.
     */
    BasicBuilder& set_top(int t = Defaults::top) {
        top = t;
        return *this;
    }

    /**
     * @param a Whether to use an approximate method to quickly find the quantile, see `Classifier::set_approximate()`.
     *
     * @return A reference to this `BasicBuilder` object.
     */
    BasicBuilder& set_approximate(bool a = Defaults::approximate) {
        approximate = a;
        return *this;
    }

    /**
     * @param n Number of threads to use.
     *
     * @return A reference to this `BasicBuilder` object.
     */
    BasicBuilder& set_num_threads(int n = Defaults::num_threads) {
        nthreads = n;
        return *this;
    }

private:
    std::vector<Reference> build_internal(const tatami::Matrix<double, int>* ref, const int* labels, const std::vector<int>& subset) const {
        std::vector<Reference> subref;
        if (approximate) {
            subref = build_indices(
                ref, 
                labels, 
                subset, 
                [](size_t nr, size_t nc, const double* ptr) { 
                    return std::shared_ptr<knncolle::Base<int, double> >(new knncolle::AnnoyEuclidean<int, double>(nr, nc, ptr)); 
                },
                nthreads
            );
        } else {
            subref = build_indices(
                ref, 
                labels, 
                subset,
                [](size_t nr, size_t nc, const double* ptr) { 
                    return std::shared_ptr<knncolle::Base<int, double> >(new knncolle::KmknnEuclidean<int, double>(nr, nc, ptr)); 
                },
                nthreads
            );
        }
        return subref;
    }

public:
    /**
     * @brief Prebuilt reference that can be directly used for annotation.
     */
    struct Prebuilt {
        /**
         * @cond
         */
        Prebuilt(Markers m, std::vector<int> s, std::vector<Reference> r) :
            markers(std::move(m)), subset(std::move(s)), references(std::move(r)) {}
        /**
         * @endcond
         */

        /**
         * A vector of vectors of ranked marker genes to be used in the classification.
         * Values are indices into the `subset` vector.
         * The set of marker genes is typically a subset of those in the input `markers` in `build()`.
         */
        Markers markers;

        /**
         * The subset of features in the test/reference datasets that were used in the classification.
         * Values are row indices into the relevant matrices.
         */
        std::vector<int> subset;

        /**
         * @return Number of labels in this reference.
         */
        size_t num_labels() const {
            return references.size();
        }

        /**
         * @return Number of profiles in this reference.
         */
        size_t num_profiles() const {
            size_t n = 0;
            for (const auto& ref : references) {
                n += ref.ranked.size();
            }
            return n;
        }

        /**
         * @cond
         */
        std::vector<Reference> references;
        /**
         * @endcond
         */
    };

    /**
     * @param ref Matrix for the reference expression profiles.
     * Rows are genes while columns are samples.
     * @param[in] labels An array of length equal to the number of columns of `ref`, containing the label for each sample.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param markers A vector of vectors of ranked marker genes for each pairwise comparison between labels, see `Markers` for more details.
     *
     * @return A `Prebuilt` instance that can be used in `run()` for annotation of a test dataset.
     */
    Prebuilt run(const tatami::Matrix<double, int>* ref, const int* labels, Markers markers) const {
        auto subset = subset_markers(markers, top);
        auto subref = build_internal(ref, labels, subset);
        return Prebuilt(std::move(markers), std::move(subset), std::move(subref));
    }

public:
    /**
     * @brief Prebuilt reference requiring an intersection of features. 
     */
    struct PrebuiltIntersection {
        /**
         * @cond
         */
        PrebuiltIntersection(Markers m, std::vector<int> mats, std::vector<int> refs, std::vector<Reference> r) :
            markers(std::move(m)), mat_subset(std::move(mats)), ref_subset(std::move(refs)), references(std::move(r)) {}
        /**
         * @endcond
         */

        /**
         * A vector of vectors of ranked marker genes to be used in the classification.
         * Values are indices into the `mat_subset` and `ref_subset` vectors for the respective matrices.
         * The set of marker genes is typically a subset of those in the input `markers` in `build()`.
         */
        Markers markers;

        /**
         * Row indices of test dataset, specifying the features in the intersection.
         * This has the same length as `ref_subset`, where corresponding entries refer to the same features in the respective datasets.
         */
        std::vector<int> mat_subset;

        /**
         * Row indices of reference dataset, specifying the features in the intersection.
         * This has the same length as `mat_subset`, where corresponding entries refer to the same features in the respective datasets.
         */
        std::vector<int> ref_subset;

        /**
         * @return Number of labels in this reference.
         */
        size_t num_labels() const {
            return references.size();
        }

        /**
         * @return Number of profiles in this reference.
         */
        size_t num_profiles() const {
            size_t n = 0;
            for (const auto& ref : references) {
                n += ref.ranked.size();
            }
            return n;
        }

        /**
         * @cond
         */
        std::vector<Reference> references;
        /**
         * @endcond
         */
    };

    /**
     * @param intersection Vector defining the intersection of features betweent the test and reference datasets.
     * Each entry is a pair where the first element is the row index in the test matrix,
     * and the second element is the row index for the corresponding feature in the reference matrix.
     * Each row index for either matrix should occur no more than once in `intersection`.
     * @param ref An expression matrix for the reference expression profiles, where rows are genes and columns are cells.
     * This should have non-zero columns.
     * @param[in] labels An array of length equal to the number of columns of `ref`, containing the label for each sample.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param markers A vector of vectors of ranked marker genes for each pairwise comparison between labels, see `Markers` for more details.
     *
     * @return A `PrebuiltIntersection` instance that can be used in `run()` for annotation of a test dataset with the same order of genes as specified in `mat_id`.
     *
     * This method deals with the case where the genes are not in the same order and number across the test and reference datasets.
     * It finds the intersection of genes and then prepares the references accordingly.
     */
    PrebuiltIntersection run(
        std::vector<std::pair<int, int> > intersection,
        const tatami::Matrix<double, int>* ref, 
        const int* labels,
        Markers markers)
    const {
        // Sorting it if it wasn't already.
        for (size_t i = 1, end = intersection.size(); i < end; ++i) {
            if (intersection[i] < intersection[i-1]) {
                std::sort(intersection.begin(), intersection.end());
                break;
            }
        }

        subset_markers(intersection, markers, top);
        auto pairs = unzip(intersection);
        auto subref = build_internal(ref, labels, pairs.second);
        return PrebuiltIntersection(std::move(markers), std::move(pairs.first), std::move(pairs.second), std::move(subref));
    }

    /**
     * @tparam Id Gene identifier for each row.
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
     * @param markers A vector of vectors of ranked marker genes for each pairwise comparison between labels, see `Markers` for more details.
     *
     * @return A `PrebuiltIntersection` instance that can be used in `run()` for annotation of a test dataset with the same order of genes as specified in `mat_id`.
     *
     * This method deals with the case where the genes are not in the same order and number across the test and reference datasets.
     * It finds the intersection of genes and then prepares the references accordingly.
     */
    template<class Id>
    PrebuiltIntersection run(
        size_t mat_nrow,
        const Id* mat_id, 
        const tatami::Matrix<double, int>* ref, 
        const Id* ref_id, 
        const int* labels,
        Markers markers) 
    const {
        auto intersection = intersect_features(mat_nrow, mat_id, ref->nrow(), ref_id);
        return run(std::move(intersection), ref, labels, std::move(markers));
    }
};

}

#endif
