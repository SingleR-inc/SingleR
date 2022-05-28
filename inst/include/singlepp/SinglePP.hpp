#ifndef SINGLEPP_SINGLEPP_HPP
#define SINGLEPP_SINGLEPP_HPP

#include "knncolle/knncolle.hpp"
#include "build_indices.hpp"
#include "annotate_cells.hpp"
#include "process_features.hpp"

#include <vector> 
#include <stdexcept>

/**
 * @file SinglePP.hpp
 *
 * @brief Defines the `SinglePP` class.
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
class SinglePP {
public:
    /**
     * @brief Default parameters for annotation.
     */
    struct Defaults {
        /**
         * See `set_quantile()` for details.
         */
        static constexpr double quantile = 0.8;

        /**
         * See `set_fine_tune_threshold()` for details.
         */
        static constexpr double fine_tune_threshold = 0.05;

        /**
         * See `set_fine_tune()` for details.
         */
        static constexpr bool fine_tune = true;

        /**
         * See `set_top()` for details.
         */
        static constexpr int top = 20;

        /**
         * See `set_approximate()` for details.
         */
        static constexpr bool approximate = false;
    };

private:
    double quantile = Defaults::quantile;
    double fine_tune_threshold = Defaults::fine_tune_threshold;
    bool fine_tune = Defaults::fine_tune;
    int top = Defaults::top;
    bool approximate = Defaults::approximate;

public:
    /**
     * @param q Quantile to use to compute a per-label score from the correlations.
     *
     * @return A reference to this `SinglePP` object.
     *
     * Values of `q` closer to 0.5 focus on the behavior of the majority of a label's reference profiles.
     * Smaller values will be more sensitive to the presence of a subset of profiles that are more similar to the test cell,
     * which can be useful when the reference profiles themselves are heterogeneous.
     */
    SinglePP& set_quantile(double q = Defaults::quantile) {
        quantile = q;
        return *this;
    }

    /**
     * @param t Threshold to use to select the top-scoring subset of labels during fine-tuning.
     * Larger values increase the chance of recovering the correct label at the cost of computational time.
     *
     * @return A reference to this `SinglePP` object.
     *
     * Needless to say, one should not set `t` to a value that is too large.
     * Otherwise, the first fine-tuning iteration would just contain all labels and there would be no reduction of the marker space.
     */
    SinglePP& set_fine_tune_threshold(double t = Defaults::fine_tune_threshold) {
        fine_tune_threshold = t;
        return *this;
    }

    /**
     * @param f Whether to perform fine-tuning.
     * This can be disabled to improve speed at the cost of accuracy.
     *
     * @return A reference to this `SinglePP` object.
     */
    SinglePP& set_fine_tune(bool f = Defaults::fine_tune) {
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
     * @return A reference to this `SinglePP` object.
     */
    SinglePP& set_top(int t = Defaults::top) {
        top = t;
        return *this;
    }

    /**
     * @param a Whether to use an approximate method to quickly find the quantile.
     * This sacrifices some accuracy for speed when labels have many reference profiles.
     *
     * @return A reference to this `SinglePP` object.
     */
    SinglePP& set_approximate(bool a = Defaults::approximate) {
        approximate = a;
        return *this;
    }

private:
    std::vector<Reference> build_internal(const tatami::Matrix<double, int>* ref, const int* labels, const std::vector<int>& subset) {
        std::vector<Reference> subref;
        if (approximate) {
            subref = build_indices(ref, labels, subset, 
                [](size_t nr, size_t nc, const double* ptr) { 
                    return std::shared_ptr<knncolle::Base<int, double> >(new knncolle::AnnoyEuclidean<int, double>(nr, nc, ptr)); 
                }
            );
        } else {
            subref = build_indices(ref, labels, subset,
                [](size_t nr, size_t nc, const double* ptr) { 
                    return std::shared_ptr<knncolle::Base<int, double> >(new knncolle::KmknnEuclidean<int, double>(nr, nc, ptr)); 
                }
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
    Prebuilt build(const tatami::Matrix<double, int>* ref, const int* labels, Markers markers) {
        auto subset = subset_markers(markers, top);
        auto subref = build_internal(ref, labels, subset);
        return Prebuilt(std::move(markers), std::move(subset), std::move(subref));
    }

public:
    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * This may have a different ordering of genes compared to the reference matrix used to create `built`,
     * provided that all genes corresponding to `built.subset` are present.
     * @param built An object produced by `build()`.
     * @param[in] mat_subset Pointer to an array of length equal to that of `built.subset`,
     * containing the index of the row of `mat` corresponding to each gene in `built.subset`.
     * That is, row `mat_subset[i]` in `mat` should be the same gene as row `built.subset[i]` in the reference matrix.
     * @param[out] best Pointer to an array of length equal to the number of columns in `mat`.
     * This is filled with the index of the assigned label for each cell.
     * @param[out] scores Vector of pointers to arrays of length equal to the number of columns in `mat`.
     * This is filled with the (non-fine-tuned) score for each label for each cell.
     * Any pointer may be `NULL` in which case the scores for that label will not be saved.
     * @param[out] delta Pointer to an array of length equal to the number of columns in `mat`.
     * This is filled with the difference between the highest and second-highest scores, possibly after fine-tuning.
     * This may also be `NULL` in which case the deltas are not reported.
     * @param already_subset Whether the rows of `mat` have already been subsetted to match `built.subset`.
     *
     * @return `best`, `scores` and `delta` are filled with their output values.
     */
    void run(const tatami::Matrix<double, int>* mat, const Prebuilt& built, const int* mat_subset, int* best, std::vector<double*>& scores, double* delta) {
        annotate_cells_simple(
            mat, 
            built.subset.size(), 
            mat_subset, 
            built.references, 
            built.markers, 
            quantile, 
            fine_tune, 
            fine_tune_threshold, 
            best, 
            scores, 
            delta
        );
        return;
    }

    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * This should have the same order and identity of genes as the reference matrix used to create `built`.
     * @param built An object produced by `build()`.
     * @param[out] best Pointer to an array of length equal to the number of columns in `mat`.
     * This is filled with the index of the assigned label for each cell.
     * @param[out] scores Vector of pointers to arrays of length equal to the number of columns in `mat`.
     * This is filled with the (non-fine-tuned) score for each label for each cell.
     * Any pointer may be `NULL` in which case the scores for that label will not be saved.
     * @param[out] delta Pointer to an array of length equal to the number of columns in `mat`.
     * This is filled with the difference between the highest and second-highest scores, possibly after fine-tuning.
     * This may also be `NULL` in which case the deltas are not reported.
     * @param already_subset Whether the rows of `mat` have already been subsetted to match `built.subset`.
     *
     * @return `best`, `scores` and `delta` are filled with their output values.
     */
    void run(const tatami::Matrix<double, int>* mat, const Prebuilt& built, int* best, std::vector<double*>& scores, double* delta) {
        run(mat, built, built.subset.data(), best, scores, delta);
        return;
    }

    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * @param ref An expression matrix for the reference expression profiles.
     * This should have non-zero columns and the same number of rows (i.e., genes) at `mat`.
     * @param[in] labels An array of length equal to the number of columns of `ref`, containing the label for each sample.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param markers A vector of vectors of ranked marker genes for each pairwise comparison between labels, see `Markers` for more details.
     * @param[out] best Pointer to an array of length equal to the number of columns in `mat`.
     * This is filled with the index of the assigned label for each cell.
     * @param[out] scores Vector of pointers of length equal to the number of labels.
     * Each pointer should point to an array of length equal to the number of columns in `mat`.
     * This is filled with the (non-fine-tuned) score for that label for each cell.
     * Any pointer may be `NULL` in which case the scores for that label will not be saved.
     * @param[out] delta Pointer to an array of length equal to the number of columns in `mat`.
     * This is filled with the difference between the highest and second-highest scores, possibly after fine-tuning.
     * This may also be `NULL` in which case the deltas are not reported.
     *
     * @return `best`, `scores` and `delta` are filled with their output values.
     */
    void run(const tatami::Matrix<double, int>* mat, const tatami::Matrix<double, int>* ref, const int* labels, Markers markers, int* best, std::vector<double*>& scores, double* delta) {
        auto prebuilt = build(ref, labels, std::move(markers));
        run(mat, prebuilt, best, scores, delta);
        return;
    }

public:
    /**
     * @brief Results of the automated annotation.
     */
    struct Results {
        /**
         * @cond
         */
        Results(size_t ncells, size_t nlabels) : best(ncells), scores(nlabels, std::vector<double>(ncells)), delta(ncells) {}

        std::vector<double*> scores_to_pointers() {
            std::vector<double*> output(scores.size());
            for (size_t s = 0; s < scores.size(); ++s) {
                output[s] = scores[s].data();
            }
            return output;
        };
        /**
         * @endcond
         */

        /** 
         * Vector of length equal to the number of cells in the test dataset,
         * containing the index of the assigned label for each cell.
         */
        std::vector<int> best;

        /**
         * Vector of length equal to the number of labels,
         * containing vectors of length equal to the number of cells in the test dataset.
         * Each vector corresponds to a label and contains the (non-fine-tuned) score for each cell.
         */
        std::vector<std::vector<double> > scores;


        /** 
         * Vector of length equal to the number of cells in the test dataset.
         * This contains the difference between the highest and second-highest scores for each cell, possibly after fine-tuning.
         */
        std::vector<double> delta;
    };

    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * If `already_subset = false`, this should have the same identity of genes as the reference matrices used in `build()`;
     * otherwise, each row should correspond to an element of `built.subset`.
     * @param built An object produced by `build()`.
     *
     * @return A `Results` object containing the assigned labels and scores.
     */
    Results run(const tatami::Matrix<double, int>* mat, const Prebuilt& built) {
        size_t nlabels = built.references.size();
        Results output(mat->ncol(), nlabels);
        auto scores = output.scores_to_pointers();
        run(mat, built, output.best.data(), scores, output.delta.data());
        return output;
    }

    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * This may have a different ordering of genes compared to the reference matrix used in `build()`,
     * provided that all genes corresponding to `built.subset` are present.
     * @param built An object produced by `build()`.
     * @param[in] mat_subset Pointer to an array of length equal to that of `built.subset`,
     * containing the index of the row of `mat` corresponding to each gene in `built.subset`.
     *
     * @return A `Results` object containing the assigned labels and scores.
     */
    Results run(const tatami::Matrix<double, int>* mat, const Prebuilt& built, const int* mat_subset) {
        size_t nlabels = built.references.size();
        Results output(mat->ncol(), nlabels);
        auto scores = output.scores_to_pointers();
        run(mat, built, mat_subset, output.best.data(), scores, output.delta.data());
        return output;
    }

    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * @param ref An expression matrix for the reference expression profiles.
     * This should have non-zero columns and the same number of rows (i.e., genes) at `mat`.
     * @param[in] labels An array of length equal to the number of columns of `ref`, containing the label for each sample.
     * The smallest label should be 0 and the largest label should be equal to the total number of unique labels minus 1.
     * @param markers A vector of vectors of ranked marker genes for each pairwise comparison between labels, see `Markers` for more details.
     *
     * @return A `Results` object containing the assigned labels and scores.
     */
    Results run(const tatami::Matrix<double, int>* mat, const tatami::Matrix<double, int>* ref, const int* labels, Markers markers) {
        auto prebuilt = build(ref, labels, std::move(markers));
        return run(mat, prebuilt);
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
         * This has the same length as `ref_indices`.
         */
        std::vector<int> mat_subset;

        /**
         * Row indices of reference dataset, specifying the features in the intersection.
         * This has the same length as `mat_indices`.
         */
        std::vector<int> ref_subset;

        /**
         * @return Number of labels in this reference.
         */
        size_t num_labels() const {
            return references.size();
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
     * @tparam Id Gene identifier for each row.
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
     * @param markers A vector of vectors of ranked marker genes for each pairwise comparison between labels, see `Markers` for more details.
     *
     * @return A `PrebuiltIntersection` instance that can be used in `run()` for annotation of a test dataset with the same order of genes as specified in `mat_id`.
     *
     * This method deals with the case where the genes are not in the same order and number across the test and reference datasets.
     * It finds the intersection of genes and then prepares the references accordingly.
     */
    template<class Id>
    PrebuiltIntersection build(
        size_t mat_nrow,
        const Id* mat_id, 
        const tatami::Matrix<double, int>* ref, 
        const Id* ref_id, 
        const int* labels,
        Markers markers) 
    {
        auto intersection = intersect_features(mat_nrow, mat_id, ref->nrow(), ref_id);
        subset_markers(intersection, markers, top);
        auto pairs = unzip(intersection);
        auto subref = build_internal(ref, labels, pairs.second);
        return PrebuiltIntersection(std::move(markers), std::move(pairs.first), std::move(pairs.second), std::move(subref));
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
     * This is filled with the index of the assigned label for each cell.
     * @param[out] scores Vector of pointers of length equal to the number of labels.
     * Each pointer should point to an array of length equal to the number of columns in `mat`.
     * This is filled with the (non-fine-tuned) score for that label for each cell.
     * Any pointer may be `NULL` in which case the scores for that label will not be saved.
     * @param[out] delta Pointer to an array of length equal to the number of columns in `mat`.
     * This is filled with the difference between the highest and second-highest scores, possibly after fine-tuning.
     * This may also be `NULL` in which case the deltas are not reported.
     * 
     * @return `best`, `scores` and `delta` are filled with their output values.
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
    {
        auto built = build(mat->nrow(), mat_id, ref, ref_id, labels, std::move(markers));
        run(mat, built, best, scores, delta);
        return;
    }

    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * @param built An object produced by `build()` with intersections.
     * @param[out] best Pointer to an array of length equal to the number of columns in `mat`.
     * This is filled with the index of the assigned label for each cell.
     * @param[out] scores Vector of pointers to arrays of length equal to the number of columns in `mat`.
     * This is filled with the (non-fine-tuned) score for each label for each cell.
     * Any pointer may be `NULL` in which case the scores for that label will not be saved.
     * @param[out] delta Pointer to an array of length equal to the number of columns in `mat`.
     * This is filled with the difference between the highest and second-highest scores, possibly after fine-tuning.
     * This may also be `NULL` in which case the deltas are not reported.
     * 
     * @return `best`, `scores` and `delta` are filled with their output values.
     */
    void run(
        const tatami::Matrix<double, int>* mat, 
        const PrebuiltIntersection& built,
        int* best,
        std::vector<double*>& scores,
        double* delta) 
    {
        annotate_cells_simple(mat, 
            built.mat_subset.size(), 
            built.mat_subset.data(), 
            built.references, 
            built.markers, 
            quantile, 
            fine_tune, 
            fine_tune_threshold, 
            best, 
            scores, 
            delta
        );
        return;
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
     *
     * @return A `Results` object containing the assigned labels and scores.
     */ 
    template<class Id>
    Results run(const tatami::Matrix<double, int>* mat, const Id* mat_id, const tatami::Matrix<double, int>* ref, const Id* ref_id, const int* labels, Markers markers) {
        auto built = build(mat->nrow(), mat_id, ref, ref_id, labels, std::move(markers));

        size_t nlabels = built.references.size();
        Results output(mat->ncol(), nlabels);
        auto scores = output.scores_to_pointers();

        run(mat, built, output.best.data(), scores, output.delta.data());
        return output;
    }

    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * @param built An object produced by `build()` with intersections.
     *
     * @return A `Results` object containing the assigned labels and scores.
     */ 
    Results run(const tatami::Matrix<double, int>* mat, const PrebuiltIntersection& built) {
        size_t nlabels = built.references.size();
        Results output(mat->ncol(), nlabels);
        auto scores = output.scores_to_pointers();

        run(mat, built, output.best.data(), scores, output.delta.data());
        return output;
    }
};

}

#endif
