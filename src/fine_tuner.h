#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "compute_scores.h"

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <tuple>

typedef std::tuple<int, double, double> tuned_stats;

class fine_tuner {
public:
    fine_tuner(size_t ngenes) : holder_left(ngenes), holder_right(ngenes), 
        scaled_left(ngenes), scaled_right(ngenes), collected(ngenes) {}

    template<class PICKER> 
    tuned_stats assign(int i, beachmat::numeric_matrix* exprs, Rcpp::NumericMatrix scores,
        const matrix_list& references, double quantile, double tune_thresh, const PICKER& commonFUN) 
    {
        exprs->get_col(i, holder_left.begin());
        auto cur_scores=scores.column(i);
        if (cur_scores.size()==0) {
            return tuned_stats(NA_INTEGER, R_NaReal, R_NaReal);
        }

        size_t topI=std::max_element(cur_scores.begin(), cur_scores.end()) - cur_scores.begin();
        int best_label=topI;
        double max_score=cur_scores[topI];
        if (cur_scores.size()==1) {
            return tuned_stats(best_label, max_score, R_NaReal);
        }

        constexpr double DUMMY=-1000;
        double next_score=DUMMY;

        labels_in_use.clear();
        double threshold=max_score - tune_thresh;
        for (size_t i=0; i<cur_scores.size(); ++i) {
            const auto& val=cur_scores[i];
            if (val >= threshold) {
                labels_in_use.push_back(i);
            }
            if (i!=topI && next_score < val) {
                next_score=val;
            }
        }

        // Check if it's unchanged, to avoid an infinite loop if the correlations are still close after fine tuning.
        bool unchanged=false;
        while (labels_in_use.size() > 1 && !unchanged) {
            commonFUN(labels_in_use, genes_in_use);

            new_scores.resize(labels_in_use.size());
            compute_scores(holder_left, 
                references, 
                labels_in_use,
                genes_in_use,
                quantile,
                holder_right,
                scaled_left,
                scaled_right,
                collected,
                all_correlations,
                new_scores.data()
            );

            size_t topI=std::max_element(new_scores.begin(), new_scores.end()) - new_scores.begin();
            best_label=labels_in_use[topI];
            max_score=new_scores[topI];
            threshold=max_score - tune_thresh;
            next_score=DUMMY;

            unchanged=true;
            for (size_t i=0; i<new_scores.size(); ++i) {
                const auto& val=new_scores[i];
                if (val >= threshold) {
                    next_labels.push_back(labels_in_use[i]);
                } else {
                    unchanged=false;
                }
                if (i!=topI && next_score < val) {
                    next_score=val;
                }
            }
            labels_in_use.swap(next_labels);
            next_labels.clear();
        }

        return tuned_stats(best_label, max_score, next_score);
    }

private:
    Rcpp::NumericVector holder_left, holder_right;
    std::vector<int> labels_in_use, next_labels, genes_in_use;

    // Entities used inside 'get_scores'.
    std::vector<double> new_scores;
    std::vector<double> scaled_left, scaled_right;
    ranked_vector collected;
    std::vector<double> all_correlations;
};
