#ifndef SINGLEPP_CLASSIFIER_HPP
#define SINGLEPP_CLASSIFIER_HPP

#include "macros.hpp"

#include "tatami/tatami.hpp"

#include "BasicBuilder.hpp"
#include "BasicScorer.hpp"

#include <vector> 
#include <stdexcept>

/**
 * @file Classifier.hpp
 *
 * @brief Defines the `Classifier` class.
 */

namespace singlepp {

/**
 * @brief Automatically assign cell type labels based on an expression matrix.
 *
 * This implements the [**SingleR**](https://bioconductor.org/packages/SingleR) algorithm for automated annotation of single-cell RNA-seq data.
 * For each cell, we compute the Spearman rank correlation between that cell and the reference expression profiles.
 * This is done using only the subset of genes that are label-specific markers,
 * most typically the top genes from pairwise comparisons between each label's expression profiles.
 * For each label, we take the correlations involving that label's reference profiles and convert it into a score.
 * The label with the highest score is used as an initial label for that cell.
 *
 * For each cell, we apply fine-tuning iterations to improve the label accuracy by refining the feature space.
 * At each iteration, we find the subset of labels with scores that are close to the maximum score according to some threshold.
 * We recompute the scores based on the markers for this label subset, and we repeat the process until only one label is left in the subset or the subset is unchanged.
 * At the end of the iterations, the label with the highest score (or the only label, if just one is left) is used as the label for the cell.
 * This process aims to remove noise by eliminating irrelevant genes when attempting to distinguish closely related labels.
 * 
 * Each label's score is defined as a user-specified quantile of the distribution of correlations across all reference profiles assigned to that label.
 * (We typically consider a large quantile, e.g., the 80% percentile of the correlations.)
 * The use of a quantile avoids problems with differences in the number of reference profiles per label;
 * in contrast, just using the "top X correlations" would implicitly favor labels with more reference profiles.
 *
 * The choice of Spearman's correlation provides some robustness against batch effects when comparing reference and test datasets.
 * Only the relative expression _within_ each cell needs to be comparable, not their relative expression across cells.
 * As a result, it does not matter whether raw counts are supplied or log-transformed expression values, as the latter is a monotonic transformation of the latter (within each cell).
 * The algorithm is also robust to differences in technologies between reference and test profiles, though it is preferable to have like-for-like comparisons. 
 *
 * @see
 * Aran D et al. (2019). 
 * Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage.
 * _Nat. Immunol._ 20, 163-172
 */
class Classifier {
public:
    /**
     * @brief Default parameters for annotation.
     */
    struct Defaults {
        /**
         * See `set_quantile()` for details.
         */
        static constexpr double quantile = BasicScorer::Defaults::quantile;

        /**
         * See `set_fine_tune_threshold()` for details.
         */
        static constexpr double fine_tune_threshold = BasicScorer::Defaults::fine_tune_threshold;

        /**
         * See `set_fine_tune()` for details.
         */
        static constexpr bool fine_tune = BasicScorer::Defaults::fine_tune;

        /**
         * See `set_top()` for details.
         */
        static constexpr int top = BasicBuilder::Defaults::top;

        /**
         * See `set_approximate()` for details.
         */
        static constexpr bool approximate = BasicBuilder::Defaults::approximate;

        /**
         * See `set_num_threads()` for details.
         */
        static constexpr int num_threads = 1;
    };

private:
    int top = Defaults::top;
    bool approximate = Defaults::approximate;
    double quantile = Defaults::quantile;
    double fine_tune_threshold = Defaults::fine_tune_threshold;
    bool fine_tune = Defaults::fine_tune;
    int nthreads = Defaults::num_threads;

public:
    /**
     * @param q Quantile to use to compute a per-label score from the correlations.
     *
     * @return A reference to this `Classifier` object.
     *
     * Values of `q` closer to 0.5 focus on the behavior of the majority of a label's reference profiles.
     * Smaller values will be more sensitive to the presence of a subset of profiles that are more similar to the test cell,
     * which can be useful when the reference profiles themselves are heterogeneous.
     */
    Classifier& set_quantile(double q = Defaults::quantile) {
        quantile = q;
        return *this;
    }

    /**
     * @param t Threshold to use to select the top-scoring subset of labels during fine-tuning.
     * Larger values increase the chance of recovering the correct label at the cost of computational time.
     *
     * @return A reference to this `Classifier` object.
     *
     * Needless to say, one should not set `t` to a value that is too large.
     * Otherwise, the first fine-tuning iteration would just contain all labels and there would be no reduction of the marker space.
     */
    Classifier& set_fine_tune_threshold(double t = Defaults::fine_tune_threshold) {
        fine_tune_threshold = t;
        return *this;
    }

    /**
     * @param f Whether to perform fine-tuning.
     * This can be disabled to improve speed at the cost of accuracy.
     *
     * @return A reference to this `Classifier` object.
     */
    Classifier& set_fine_tune(bool f = Defaults::fine_tune) {
        fine_tune = f;
        return *this;
    }

    /**
     * @param t Number of top markers to use from each pairwise comparison between labels.
     * Larger values improve the stability of the correlations at the cost of increasing noise and computational work.
     *
     * Setting it to a negative value will instruct `run()` to use all supplied markers.
     * This is useful in situations where the supplied markers have already been curated.
     *
     * @return A reference to this `Classifier` object.
     */
    Classifier& set_top(int t = Defaults::top) {
        top = t;
        return *this;
    }

    /**
     * @param a Whether to use an approximate method to quickly find the quantile.
     * This sacrifices some accuracy for speed when labels have many reference profiles.
     *
     * @return A reference to this `Classifier` object.
     */
    Classifier& set_approximate(bool a = Defaults::approximate) {
        approximate = a;
        return *this;
    }

    /**
     * @param n Number of threads to use.
     *
     * @return A reference to this `Classifier` object.
     */
    Classifier& set_num_threads(int n = Defaults::num_threads) {
        nthreads = n;
        return *this;
    }

private:
    BasicBuilder::Prebuilt build_reference(const tatami::Matrix<double, int>* ref, const int* labels, Markers markers) const {
        BasicBuilder builder;
        builder
            .set_top(top)
            .set_approximate(approximate)
            .set_num_threads(nthreads);
        return builder.run(ref, labels, std::move(markers));
    }

    template<class Id>
    BasicBuilder::PrebuiltIntersection build_reference(size_t mat_nrow, const Id* mat_id, const tatami::Matrix<double, int>* ref, const Id* ref_id, const int* labels, Markers markers) const {
        BasicBuilder builder;
        builder
            .set_top(top)
            .set_approximate(approximate)
            .set_num_threads(nthreads);
        return builder.run(mat_nrow, mat_id, ref, ref_id, labels, std::move(markers));
    }

    BasicScorer set_up_scorer() const {
        BasicScorer scorer;
        scorer
            .set_quantile(quantile)
            .set_fine_tune(fine_tune)
            .set_fine_tune_threshold(fine_tune_threshold)
            .set_num_threads(nthreads);
        return scorer;
    }

public:
    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * @param ref An expression matrix for the reference expression profiles.
     * This should have non-zero columns and the same number of rows (i.e., genes) at `mat`.
     * @param[in] labels An array of length equal to the number of columns of `ref`, containing the label for each sample.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param markers A vector of vectors of ranked marker genes for each pairwise comparison between labels, see `Markers` for more details.
     * @param[out] best Pointer to an array of length equal to the number of columns in `mat`.
     * On output, this is filled with the index of the assigned label for each cell.
     * @param[out] scores Vector of pointers of length equal to the number of labels.
     * Each pointer should point to an array of length equal to the number of columns in `mat`.
     * On output, this is filled with the (non-fine-tuned) score for that label for each cell.
     * Any pointer may be `NULL` in which case the scores for that label will not be reported.
     * @param[out] delta Pointer to an array of length equal to the number of columns in `mat`.
     * On output, this is filled with the difference between the highest and second-highest scores, possibly after fine-tuning.
     * This may also be `NULL` in which case the deltas are not reported.
     */
    void run(const tatami::Matrix<double, int>* mat, const tatami::Matrix<double, int>* ref, const int* labels, Markers markers, int* best, std::vector<double*>& scores, double* delta) const {
        auto prebuilt = build_reference(ref, labels, std::move(markers));
        set_up_scorer().run(mat, prebuilt, best, scores, delta);
        return;
    }

    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * @param ref An expression matrix for the reference expression profiles.
     * This should have non-zero columns and the same number of rows (i.e., genes) at `mat`.
     * @param[in] labels An array of length equal to the number of columns of `ref`, containing the label for each sample.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param markers A vector of vectors of ranked marker genes for each pairwise comparison between labels, see `Markers` for more details.
     *
     * @return A `BasicScorer::Results` object containing the assigned labels and scores.
     */
    BasicScorer::Results run(const tatami::Matrix<double, int>* mat, const tatami::Matrix<double, int>* ref, const int* labels, Markers markers) const {
        auto prebuilt = build_reference(ref, labels, std::move(markers));
        return set_up_scorer().run(mat, prebuilt);
    }

public:
    /**
     * @tparam Id Gene identifier for each row.
     *
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * @param[in] mat_id Pointer to an array of identifiers of length equal to the number of rows of `mat`.
     * This should contain a unique identifier for each row of `mat` (typically a gene name or index).
     * @param ref An expression matrix for the reference expression profiles, where rows are genes and columns are cells.
     * This should have non-zero columns.
     * @param[in] ref_id Pointer to an array of identifiers of length equal to the number of rows of any `ref`.
     * This should contain a unique identifier for each row in `ref`, and should be comparable to `mat_id`.
     * @param[in] labels An array of length equal to the number of columns of `ref`, containing the label for each sample.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param markers A vector of vectors of ranked marker genes for each pairwise comparison between labels, see `Markers` for more details.
     * @param[out] best Pointer to an array of length equal to the number of columns in `mat`.
     * On output, this is filled with the index of the assigned label for each cell.
     * @param[out] scores Vector of pointers of length equal to the number of labels.
     * Each pointer should point to an array of length equal to the number of columns in `mat`.
     * On output, this is filled with the (non-fine-tuned) score for that label for each cell.
     * Any pointer may be `NULL` in which case the scores for that label will not be reported.
     * @param[out] delta Pointer to an array of length equal to the number of columns in `mat`.
     * On output, this is filled with the difference between the highest and second-highest scores, possibly after fine-tuning.
     * This may also be `NULL` in which case the deltas are not reported.
     * 
     * This version of `run()` applies an intersection to find the common genes between `mat` and `ref`, based on their shared values in `mat_id` and `ref_id`.
     * The annotation is then performed using only the subset of common genes.
     * The aim is to easily accommodate differences in feature annotation between the test and reference profiles.
     */
    template<class Id>
    void run(
        const tatami::Matrix<double, int>* mat, 
        const Id* mat_id, 
        const tatami::Matrix<double, int>* ref, 
        const Id* ref_id, 
        const int* labels,
        Markers markers, 
        int* best,
        std::vector<double*>& scores,
        double* delta) 
    const {
        auto built = build_reference(mat->nrow(), mat_id, ref, ref_id, labels, std::move(markers));
        set_up_scorer().run(mat, built, best, scores, delta);
        return;
    }

    /**
     * @tparam Id Gene identifier for each row.
     *
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * @param[in] mat_id Pointer to an array of identifiers of length equal to the number of rows of `mat`.
     * This should contain a unique identifier for each row of `mat` (typically a gene name or index).
     * @param ref An expression matrix for the reference expression profiles, where rows are genes and columns are cells.
     * This should have non-zero columns.
     * @param[in] ref_id Pointer to an array of identifiers of length equal to the number of rows of any `ref`.
     * This should contain a unique identifier for each row in `ref`, and should be comparable to `mat_id`.
     * @param[in] labels An array of length equal to the number of columns of `ref`, containing the label for each sample.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param markers A vector of vectors of ranked marker genes for each pairwise comparison between labels, see `Markers` for more details.
     *
     * @return A `BasicScorer::Results` object containing the assigned labels and scores.
     */ 
    template<class Id>
    BasicScorer::Results run(const tatami::Matrix<double, int>* mat, const Id* mat_id, const tatami::Matrix<double, int>* ref, const Id* ref_id, const int* labels, Markers markers) const {
        auto built = build_reference(mat->nrow(), mat_id, ref, ref_id, labels, std::move(markers));
        return set_up_scorer().run(mat, built);
    }
};

}

#endif
