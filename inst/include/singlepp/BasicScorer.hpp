#ifndef SINGLEPP_BASIC_SCORER_HPP
#define SINGLEPP_BASIC_SCORER_HPP

#include "macros.hpp"

#include "tatami/tatami.hpp"

#include "annotate_cells.hpp"
#include "process_features.hpp"
#include "BasicBuilder.hpp"

#include <vector> 
#include <stdexcept>

/**
 * @file BasicScorer.hpp
 *
 * @brief Defines the `BasicScorer` class.
 */

namespace singlepp {

/**
 * @brief Classify cells from a single pre-built reference dataset.
 *
 * This class uses the pre-built reference from `BasicBuilder` to classify each cell in a test dataset.
 * The algorithm and parameters are the same as described for the `Classifier` class;
 * in fact, `Classifier` just calls `BasicBuilder::run()` and then `BasicScorer::run()`.
 *
 * It is occasionally useful to call these two functions separately if the same reference dataset is to be used multiple times,
 * e.g., on different test datasets or with different parameters.
 * In such cases, we can save time by avoiding redundant builds;
 * we just have to call `BasicScorer::run()` in all subsequent uses of the pre-built reference.
 */
class BasicScorer {
public:
    /**
     * @brief Default parameters for classification.
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
         * See `set_num_threads()` for details.
         */
        static constexpr int num_threads = 1;
    };

private:
    double quantile = Defaults::quantile;
    double fine_tune_threshold = Defaults::fine_tune_threshold;
    bool fine_tune = Defaults::fine_tune;
    int nthreads = Defaults::num_threads;

public:
    /**
     * @param q Quantile to use to compute a per-label score from the correlations, see `Classifier::set_quantile()`.
     *
     * @return A reference to this `BasicScorer` object.
     */
    BasicScorer& set_quantile(double q = Defaults::quantile) {
        quantile = q;
        return *this;
    }

    /**
     * @param t Threshold to use to select the top-scoring subset of labels during fine-tuning, see `Classifier::set_fine_tune_threshold()`.
     *
     * @return A reference to this `BasicScorer` object.
     */
    BasicScorer& set_fine_tune_threshold(double t = Defaults::fine_tune_threshold) {
        fine_tune_threshold = t;
        return *this;
    }

    /**
     * @param f Whether to perform fine-tuning, see `Classifier::set_fine_tune()`.
     *
     * @return A reference to this `BasicScorer` object.
     */
    BasicScorer& set_fine_tune(bool f = Defaults::fine_tune) {
        fine_tune = f;
        return *this;
    }

    /**
     * @param n Number of threads to use.
     *
     * @return A reference to this `BasicScorer` object.
     */
    BasicScorer& set_num_threads(int n = Defaults::num_threads) {
        nthreads = n;
        return *this;
    }

public:
    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * This may have a different ordering of genes compared to the reference matrix used to create `built`,
     * provided that all genes corresponding to `built.subset` are present.
     * @param built An object produced by `BasicBuilder::build()`.
     * @param[in] mat_subset Pointer to an array of length equal to that of `built.subset`,
     * containing the index of the row of `mat` corresponding to each gene in `built.subset`.
     * That is, row `mat_subset[i]` in `mat` should be the same gene as row `built.subset[i]` in the reference matrix.
     * @param[out] best Pointer to an array of length equal to the number of columns in `mat`.
     * On output, this is filled with the index of the assigned label for each cell.
     * @param[out] scores Vector of pointers to arrays of length equal to the number of columns in `mat`.
     * On output, this is filled with the (non-fine-tuned) score for each label for each cell.
     * Any pointer may be `NULL` in which case the scores for that label will not be reported.
     * @param[out] delta Pointer to an array of length equal to the number of columns in `mat`.
     * On output, this is filled with the difference between the highest and second-highest scores, possibly after fine-tuning.
     * This may also be `NULL` in which case the deltas are not reported.
     */
    void run(const tatami::Matrix<double, int>* mat, const BasicBuilder::Prebuilt& built, const int* mat_subset, int* best, std::vector<double*>& scores, double* delta) const {
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
            delta,
            nthreads
        );
        return;
    }

    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * This should have the same order and identity of genes as the reference matrix used to create `built`.
     * @param built An object produced by `BasicBuilder::build()`.
     * @param[out] best Pointer to an array of length equal to the number of columns in `mat`.
     * On output, this is filled with the index of the assigned label for each cell.
     * @param[out] scores Vector of pointers to arrays of length equal to the number of columns in `mat`.
     * On output, this is filled with the (non-fine-tuned) score for each label for each cell.
     * Any pointer may be `NULL` in which case the scores for that label will not be reported.
     * @param[out] delta Pointer to an array of length equal to the number of columns in `mat`.
     * On output, this is filled with the difference between the highest and second-highest scores, possibly after fine-tuning.
     * This may also be `NULL` in which case the deltas are not reported.
     */
    void run(const tatami::Matrix<double, int>* mat, const BasicBuilder::Prebuilt& built, int* best, std::vector<double*>& scores, double* delta) const {
        run(mat, built, built.subset.data(), best, scores, delta);
        return;
    }

public:
    /**
     * @brief Automated classification results.
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
     * Each row should correspond to an element of `Prebuilt::subset`.
     * @param built An object produced by `BasicBuilder::build()`.
     *
     * @return A `Results` object containing the assigned labels and scores.
     */
    Results run(const tatami::Matrix<double, int>* mat, const BasicBuilder::Prebuilt& built) const {
        size_t nlabels = built.references.size();
        Results output(mat->ncol(), nlabels);
        auto scores = output.scores_to_pointers();
        run(mat, built, output.best.data(), scores, output.delta.data());
        return output;
    }

    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * This may have a different ordering of genes compared to the reference matrix used in `build()`,
     * provided that all genes corresponding to `Prebuilt::subset` are present.
     * @param built An object produced by `BasicBuilder::build()`.
     * @param[in] mat_subset Pointer to an array of length equal to that of `Prebuilt::subset`,
     * containing the index of the row of `mat` corresponding to each gene in `Prebuilt::subset`.
     *
     * @return A `Results` object containing the assigned labels and scores.
     */
    Results run(const tatami::Matrix<double, int>* mat, const BasicBuilder::Prebuilt& built, const int* mat_subset) const {
        size_t nlabels = built.references.size();
        Results output(mat->ncol(), nlabels);
        auto scores = output.scores_to_pointers();
        run(mat, built, mat_subset, output.best.data(), scores, output.delta.data());
        return output;
    }

public:
    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * @param built An object produced by `build()` with intersections.
     * @param[out] best Pointer to an array of length equal to the number of columns in `mat`.
     * On output, this is filled with the index of the assigned label for each cell.
     * @param[out] scores Vector of pointers to arrays of length equal to the number of columns in `mat`.
     * On output, this is filled with the (non-fine-tuned) score for each label for each cell.
     * Any pointer may be `NULL` in which case the scores for that label will not be reported.
     * @param[out] delta Pointer to an array of length equal to the number of columns in `mat`.
     * On output, tkkhis is filled with the difference between the highest and second-highest scores, possibly after fine-tuning.
     * This may also be `NULL` in which case the deltas are not reported.
     */
    void run(
        const tatami::Matrix<double, int>* mat, 
        const BasicBuilder::PrebuiltIntersection& built,
        int* best,
        std::vector<double*>& scores,
        double* delta) 
    const {
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
            delta,
            nthreads
        );
        return;
    }

    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * @param built An object produced by `build()` with intersections.
     *
     * @return A `Results` object containing the assigned labels and scores.
     */ 
    Results run(const tatami::Matrix<double, int>* mat, const BasicBuilder::PrebuiltIntersection& built) const {
        size_t nlabels = built.references.size();
        Results output(mat->ncol(), nlabels);
        auto scores = output.scores_to_pointers();

        run(mat, built, output.best.data(), scores, output.delta.data());
        return output;
    }
};

}


#endif
