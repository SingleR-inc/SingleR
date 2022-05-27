#ifndef SINGLEPP_RECOMPUTE_SCORES_HPP
#define SINGLEPP_RECOMPUTE_SCORES_HPP

#include "compute_scores.hpp"
#include "scaled_ranks.hpp"
#include "IntegratedBuilder.hpp"

#include <vector>
#include <unordered_map>
#include <unordered_set>

/**
 * @file IntegratedScorer.hpp
 *
 * @brief Integrate classifications from multiple references.
 */

namespace singlepp {

/**
 * @brief Integrate classifications from multiple references.
 *
 * In situations where multiple reference datasets are available,
 * we would like to obtain a single prediction for each cell from all of those references.
 * This is somewhat tricky as the different references are likely to contain strong batch effects,
 * complicating the calculation of marker genes between labels from different references (and thus precluding direct use of the usual `SinglePP::run()`).
 * The labels themselves also tend to be inconsistent, e.g., different vocabularies and resolutions, making it difficult to define sensible groups in a combined "super-reference".
 *
 * To avoid these issues, we first perform classification within each reference individually.
 * For each test cell, we identify its predicted label from a given reference, and we collect all the marker genes for that label (across all pairwise comparisons in that reference).
 * After doing this for each reference, we pool all of the collected markers to obtain a common set of interesting genes.
 * We then compute the correlation-based score between the test cell's expression profile and its predicted label from each reference, using that common set of genes.
 * The label with the highest score is considered the best representative across all references.
 *
 * This strategy is similar to `SinglePP` without fine-tuning, except that we are choosing between the best labels from each reference rather than all labels from one reference.
 * The main idea is to create a common feature set so that the correlations can be reasonably compared across references.
 * Note that differences in the feature sets across references are tolerated by simply ignoring missing genes when computing the correlations.
 * This reduces the comparability of the scores as the effective feature set will vary a little (or a lot, depending) across references;
 * nonetheless, it is preferred to taking the intersection, which is liable to leave us with very few genes.
 *
 * Our approach avoids any direct comparison between the expression profiles of different references,
 * allowing us to side-step the question of how to deal with the batch effects.
 * Similarly, we defer responsibility on solving the issue of label heterogeneity,
 * by just passing along the existing labels and leaving it to the user's interpretation.
 */
class IntegratedScorer {
private:
    double quantile = SinglePP::Defaults::quantile;

public:
    /**
     * @param q Quantile to use to compute a per-label score from the correlations.
     *
     * @return A reference to this `SinglePP` object.
     *
     * See `SinglePP::set_quantile()` for more details.
     */
    IntegratedScorer& set_quantile(double q = SinglePP::Defaults::quantile) {
        quantile = q;
        return *this;
    }

private:
    /* Here, we've split out some of the functions for easier reading.
     * Otherwise everything would land in a single mega-function.
     */

    static void build_universe(int cell,
        const std::vector<const int*>& assigned,
        const std::vector<IntegratedReference>& references, 
        std::unordered_set<int>& uset, 
        std::vector<int>& uvec) 
    {
        uset.clear();
        for (size_t r = 0; r < references.size(); ++r) {
            auto best = assigned[r][cell];
            const auto& markers = references[r].markers[best];
            uset.insert(markers.begin(), markers.end());
        }

        uvec.clear();
        uvec.insert(uvec.end(), uset.begin(), uset.end());
        std::sort(uvec.begin(), uvec.end());
        return;
    }

    static void fill_ranks(
        const tatami::Matrix<double, int>* mat, 
        const std::vector<int>& universe, 
        int cell,
        std::vector<double>& buffer,
        tatami::Workspace* wrk,
        RankedVector<double, int>& data_ranked) 
    {
        data_ranked.clear();
        if (universe.empty()) {
            return;
        }

        size_t first = universe.front();
        size_t last = universe.back() + 1;
        auto ptr = mat->column(cell, buffer.data(), first, last, wrk);

        for (auto u : universe) {
            data_ranked.emplace_back(ptr[u - first], u);                    
        }
        std::sort(data_ranked.begin(), data_ranked.end());
        return;
    }

    void prepare_mapping(
        const IntegratedReference& ref, 
        const std::vector<int>& universe,
        std::unordered_map<int, int>& mapping)
    {
        mapping.clear();
        if (ref.check_availability) {
            // If we need to check availability, we reconstruct the mapping
            // for the intersection of features available in this reference.
            // We then calculate the scaled ranks for the data.
            int counter = 0;
            for (auto c : universe) {
                if (ref.available.find(c) != ref.available.end()) {
                    mapping[c] = counter;
                    ++counter;
                }
            }
        } else {
            // If we don't need to check availability, the mapping is
            // much simpler. Technically, we could do this in the outer
            // loop if none of the references required checks, but this
            // seems like a niche optimization in practical settings.
            for (size_t s = 0; s < universe.size(); ++s) {
                auto u = universe[s];
                mapping[u] = s;
            }
        } 
        return;
    }

public:
    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * The identity of the rows should be consistent with the arguments used in `IntegratedBuilder::add()`.
     * @param[in] assigned Vector of pointers of length equal to the number of references.
     * Each pointer should point to an array of length equal to the number of columns in `mat`,
     * containing the assigned label for each column in each reference.
     * @param references Vector of integrated references produced by `IntegratedBuilder::finish()`.
     * @param[out] best Pointer to an array of length equal to the number of columns in `mat`.
     * This is filled with the index of the reference with the best label for each cell.
     * @param[out] scores Vector of pointers of length equal to the number of references.
     * Each pointer should point to an array of length equal to the number of columns in `mat`.
     * This is filled with the (non-fine-tuned) score for the best label of that reference for each cell.
     * Any pointer may be `NULL` in which case the scores for that label will not be saved.
     * @param[out] delta Pointer to an array of length equal to the number of columns in `mat`.
     * This is filled with the difference between the highest and second-highest scores.
     * This may also be `NULL` in which case the deltas are not reported.
     *
     * @return `best`, `scores` and `delta` are filled with their output values.
     */
    void run(
        const tatami::Matrix<double, int>* mat,
        const std::vector<const int*>& assigned,
        const std::vector<IntegratedReference>& references,
        int* best,
        std::vector<double*>& scores,
        double* delta)
    {
        size_t NR = mat->nrow();
        size_t NC = mat->ncol();

        #pragma omp parallel
        {
            std::vector<double> buffer(NR);
            auto wrk = mat->new_workspace(false);

            RankedVector<double, int> data_ranked, data_ranked2;
            data_ranked.reserve(NR);
            data_ranked2.reserve(NR);
            RankedVector<int, int> ref_ranked;
            ref_ranked.reserve(NR);

            std::vector<double> scaled_data(NR);
            std::vector<double> scaled_ref(NR);
            std::unordered_set<int> universe_tmp;
            std::vector<int> universe;
            std::unordered_map<int, int> mapping;
            std::vector<double> all_correlations;

            #pragma omp for
            for (size_t i = 0; i < NC; ++i) {
                build_universe(i, assigned, references, universe_tmp, universe);
                fill_ranks(mat, universe, i, buffer, wrk.get(), data_ranked);

                // Scanning through each reference and computing the score for the best group.
                double best_score = -1000, next_best = -1000;
                int best_ref = 0;
                
                for (size_t r = 0; r < references.size(); ++r) {
                    const auto& ref = references[r];
                    prepare_mapping(ref, universe, mapping);

                    scaled_ref.resize(mapping.size());
                    scaled_data.resize(mapping.size());

                    data_ranked2.clear();
                    subset_ranks(data_ranked, data_ranked2, mapping);
                    scaled_ranks(data_ranked2, scaled_data.data());

                    // Now actually calculating the score for the best group for
                    // this cell in this reference. This assumes that 'ref.ranked'
                    // already contains sorted pairs where the indices refer to the
                    // rows of the original data matrix.
                    auto best = assigned[r][i];
                    const auto& best_ranked = ref.ranked[best];
                    all_correlations.clear();

                    for (size_t s = 0; s < best_ranked.size(); ++s) {
                        ref_ranked.clear();
                        subset_ranks(best_ranked[s], ref_ranked, mapping);
                        scaled_ranks(ref_ranked, scaled_ref.data());
                        double cor = distance_to_correlation(scaled_ref.size(), scaled_data, scaled_ref);
                        all_correlations.push_back(cor);
                    }

                    double score = correlations_to_scores(all_correlations, quantile);
                    if (scores[r]) {
                        scores[r][i] = score;
                    }
                    if (score > best_score) {
                        next_best = best_score;
                        best_score = score;
                        best_ref = r;
                    } else if (score > next_best) {
                        next_best = score;
                    }
                }

                if (best) {
                    best[i] = best_ref;
                }
                if (delta && references.size() > 1) {
                    delta[i] = best_score - next_best;
                }
            }
        }
        return;
    }

public:
    /**
     * @brief Results of the integrated annotation.
     */
    struct Results {
        /**
         * @cond
         */
        Results(size_t ncells, size_t nrefs) : best(ncells), scores(nrefs, std::vector<double>(ncells)), delta(ncells) {}

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
         * containing the index of the reference with the top-scoring label for each cell.
         */
        std::vector<int> best;

        /**
         * Vector of length equal to the number of references,
         * containing vectors of length equal to the number of cells in the test dataset.
         * Each vector corresponds to a reference and contains the score for the best label in that reference for each cell.
         */
        std::vector<std::vector<double> > scores;

        /** 
         * Vector of length equal to the number of cells in the test dataset.
         * This contains the difference between the highest and second-highest scores for each cell.
         */
        std::vector<double> delta;
    };

    /**
     * @param mat Expression matrix of the test dataset, where rows are genes and columns are cells.
     * The identity of the rows should be consistent with the arguments used in `IntegratedBuilder::add()`.
     * @param[in] assigned Vector of pointers of length equal to the number of references.
     * Each pointer should point to an array of length equal to the number of columns in `mat`,
     * containing the assigned label for each column in each reference.
     * @param references Vector of integrated references produced by `IntegratedBuilder::finish()`.
     *
     * @return A `Results` object containing the assigned labels and scores.
     */
    Results run( 
        const tatami::Matrix<double, int>* mat,
        const std::vector<const int*>& assigned,
        const std::vector<IntegratedReference>& references)
    {
        Results output(mat->ncol(), references.size());
        auto scores = output.scores_to_pointers();
        run(mat, assigned, references, output.best.data(), scores, output.delta.data());
        return output;
    }
};

}

#endif
