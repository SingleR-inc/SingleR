#ifndef SINGLEPP_ANNOTATE_CELLS_HPP
#define SINGLEPP_ANNOTATE_CELLS_HPP

#include "tatami/tatami.hpp"

#include "scaled_ranks.hpp"
#include "process_features.hpp"
#include "build_indices.hpp"
#include "fine_tune.hpp"

#include <vector>
#include <algorithm>
#include <cmath>

namespace singlepp {

inline void annotate_cells_simple(
    const tatami::Matrix<double, int>* mat,
    const std::vector<int>& subset,
    const std::vector<Reference>& ref,
    const Markers& markers,
    double quantile,
    bool fine_tune,
    double threshold,
    int* best, 
    std::vector<double*>& scores,
    double* delta) 
{
    size_t first = 0, last = 0;
    if (subset.size()) {
        // Assumes that 'subset' is sorted.
        first = subset.front();
        last = subset.back() + 1;
    }
    const size_t NC = mat->ncol();

    // Figuring out how many neighbors to keep and how to compute the quantiles.
    const size_t NL = ref.size();
    std::vector<int> search_k(NL);
    std::vector<std::pair<double, double> > coeffs(NL);
    for (size_t r = 0; r < NL; ++r) {
        double denom = ref[r].index->nobs() - 1;
        double prod = denom * (1 - quantile);
        auto k = std::ceil(prod) + 1;
        search_k[r] = k;

        // `(1 - quantile) - (k - 2) / denom` represents the gap to the smaller quantile,
        // while `(k - 1) / denom - (1 - quantile)` represents the gap from the larger quantile.
        // The size of the gap is used as the weight for the _other_ quantile, i.e., 
        // the closer you are to a quantile, the higher the weight.
        // We convert these into proportions by dividing by their sum, i.e., `1/denom`.
        coeffs[r].first = static_cast<double>(k - 1) - prod;
        coeffs[r].second = prod - static_cast<double>(k - 2);
    }

    #pragma omp parallel
    {
        std::vector<double> buffer(last - first);
        auto wrk = mat->new_workspace(false);

        RankedVector<double, int> vec(subset.size());
        std::vector<double> scaled(subset.size());

        FineTuner ft;
        std::vector<double> curscores(NL);

        #pragma omp for
        for (size_t c = 0; c < NC; ++c) {
            auto ptr = mat->column(c, buffer.data(), first, last, wrk.get());
            fill_ranks(subset, ptr, vec, first);
            scaled_ranks(vec, scaled.data());

            curscores.resize(NL);
            for (size_t r = 0; r < NL; ++r) {
                size_t k = search_k[r];
                auto current = ref[r].index->find_nearest_neighbors(scaled.data(), k);

                double last = current[k - 1].second;
                last = 1 - 2 * last * last;
                if (k == 1) {
                    curscores[r] = last;
                } else {
                    double next = current[k - 2].second;
                    next = 1 - 2 * next * next;
                    curscores[r] = coeffs[r].first * next + coeffs[r].second * last;
                }

                if (scores[r]) {
                    scores[r][c] = curscores[r];
                }
            }

            if (!fine_tune) {
                auto top = std::max_element(curscores.begin(), curscores.end());
                best[c] = top - curscores.begin();
                if (delta) {
                    if (curscores.size() > 1) {
                        double topscore = *top;
                        *top = -100;
                        delta[c] = topscore - *std::max_element(curscores.begin(), curscores.end());
                    } else {
                        delta[c] = std::numeric_limits<double>::quiet_NaN();
                    }
                }
            } else {
                auto tuned = ft.run(vec, ref, markers, curscores, quantile, threshold);
                best[c] = tuned.first;
                if (delta) {
                    delta[c] = tuned.second;
                }
            }
        }
    }

    return;
}

}

#endif
