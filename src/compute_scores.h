#ifndef COMPUTE_SCORES_H
#define COMPUTE_SCORES_H

#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "scaled_ranks.h"

#include <vector>
#include <algorithm>

typedef std::vector<std::unique_ptr<beachmat::numeric_matrix> > matrix_list;

inline void compute_scores(
    const Rcpp::NumericVector& target,
    const matrix_list& references, 
    const int* labels, 
    size_t nlabels,
    const std::vector<int>& genes, 
    double quantile,
    Rcpp::NumericVector& holding_ref,
    std::vector<double>& scaled_left,
    std::vector<double>& scaled_right,
    ranked_vector collected,
    std::vector<double>& all_correlations,
    double* new_scores
) {
    scaled_ranks(target.begin(), genes, collected, scaled_left);

    for (size_t l=0; l<nlabels; ++l) {
        auto current=references[labels[l]].get();
        const size_t ncells=current->get_ncol();
        all_correlations.clear();
        all_correlations.reserve(ncells);

        for (size_t c=0; c<ncells; ++c) {
            current->get_col(c, holding_ref.begin());
            scaled_ranks(holding_ref.begin(), genes, collected, scaled_right);

            double dist=0;
            for (size_t j=0; j<scaled_left.size(); ++j) {
                const double tmp=scaled_left[j] - scaled_right[j];
                dist+=tmp*tmp;
            }
            all_correlations.push_back(1 - 2*dist);
        }

        if (quantile==1 || ncells==1) {
            new_scores[l]=*std::max_element(all_correlations.begin(), all_correlations.end());
        } else {
            // See logic in .find_nearest_quantile().
            const double denom=ncells-1;
            const size_t qn=std::floor(denom * quantile) + 1;

            // Technically, I should do (qn-1)+1, with the first -1 being to get zero-indexed values
            // and the second +1 to obtain the ceiling. But they cancel out, so I won't.
            std::nth_element(all_correlations.begin(), all_correlations.begin()+qn, all_correlations.end());
            const double rightval=all_correlations[qn];

            // Do NOT be tempted to do the second nth_element with the end at end()+qn;
            // this does not handle ties properly.
            std::nth_element(all_correlations.begin(), all_correlations.begin()+qn-1, all_correlations.end());
            const double leftval=all_correlations[qn-1];

            const double rightweight=quantile - ((qn-1)/denom);
            const double leftweight=(qn/denom) - quantile;
            new_scores[l]=(rightval * rightweight + leftval * leftweight)/(rightweight + leftweight);
        }
    }

    return;
}

#endif
