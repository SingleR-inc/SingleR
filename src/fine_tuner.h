#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "scaled_ranks.h"

#include <set>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iostream>

typedef std::vector<std::unique_ptr<beachmat::numeric_matrix> > matrix_list;

class fine_tuner {
public:
    fine_tuner(size_t ngenes) : holder_left(ngenes), holder_right(ngenes), scaled_left(ngenes), scaled_right(ngenes), collected(ngenes) {}

    template<class PICKER> 
    int assign(int i, beachmat::numeric_matrix* exprs, Rcpp::NumericMatrix scores,
        const matrix_list& references, double quantile, double tune_thresh, PICKER commonFUN) 
    {
        exprs->get_col(i, holder_left.begin());
        auto cur_scores=scores.column(i);
        int topI=std::max_element(cur_scores.begin(), cur_scores.end()) - cur_scores.begin();
        double threshold=cur_scores[topI] - tune_thresh;

        labels_in_use.clear();
        for (size_t i=0; i<cur_scores.size(); ++i) {
            if (cur_scores[i] >= threshold) {
                labels_in_use.push_back(i);
            }
        }

        // Check if it's unchanged, to avoid an infinite loop if the correlations are still close after fine tuning.
        bool unchanged=false;
        while (labels_in_use.size() > 1 && !unchanged) {
            commonFUN(labels_in_use, genes_in_use);
            get_scores(references, quantile);
            topI=std::max_element(new_scores.begin(), new_scores.end()) - new_scores.begin();
            threshold=new_scores[topI] - tune_thresh;

            unchanged=true;
            for (size_t i=0; i<new_scores.size(); ++i) {
                if (new_scores[i] >= threshold) {
                    next_labels.push_back(labels_in_use[i]);
                } else {
                    unchanged=false;
                }
            }
            labels_in_use.swap(next_labels);
            next_labels.clear();
        }

        if (labels_in_use.size()==1L) {
            return labels_in_use.front();
        } else if (labels_in_use.size()==0L) {
            return NA_INTEGER;
        } else {
            return labels_in_use[topI];
        }
    }

    void get_scores(const matrix_list& references, double quantile) {
        scaled_ranks(holder_left.begin(), genes_in_use, collected, scaled_left);
        new_scores.clear();
        new_scores.reserve(labels_in_use.size());

        for (auto l : labels_in_use) {
            auto current=references[l].get();
            const size_t ncells=current->get_ncol();
            all_correlations.clear();
            all_correlations.reserve(ncells);

            for (size_t c=0; c<ncells; ++c) {
                current->get_col(c, holder_right.begin());
                scaled_ranks(holder_right.begin(), genes_in_use, collected, scaled_right);

                double dist=0;
                for (size_t j=0; j<scaled_left.size(); ++j) {
                    const double tmp=scaled_left[j] - scaled_right[j];
                    dist+=tmp*tmp;
                }
                all_correlations.push_back(-(1 - 2*dist)); // -1, for sorting purposes.
            }

            const int k=std::max(1, static_cast<int>((1-quantile) * ncells + 0.5)) - 1;
            std::nth_element(all_correlations.begin(), all_correlations.begin()+k, all_correlations.end());
            new_scores.push_back(-all_correlations[k]);
        }
    }
private:
    Rcpp::NumericVector holder_left, holder_right;
    std::vector<int> labels_in_use, next_labels;
    std::set<int> genes_in_use;

    // Entities used inside 'get_scores'.
    std::vector<double> new_scores;
    std::vector<double> scaled_left, scaled_right;
    ranked_vector collected;
    std::vector<double> all_correlations;
};
