#include "beachmat/numeric_matrix.h"
#include "fine_tuner.h"

struct de_markers {
    de_markers(Rcpp::List marker_genes) : collected(marker_genes.size()) {
        for (size_t i=0; i<marker_genes.size(); ++i) {
            Rcpp::List internals=marker_genes[i];
            auto& current=collected[i];
            for (size_t j=0; j<internals.size(); ++j) {
                current.push_back(Rcpp::IntegerVector(internals[j]));
            }
        }
        return;
    }
    void operator()(const std::vector<int>& labels, std::set<int>& genes) {
        genes.clear();
        for (auto l : labels) {
            for (auto l2 : labels) {
                if (l!=l2) {
                    auto& current=collected[l][l2];
                    genes.insert(current.begin(), current.end());
                }
            }
        }
        return;
    }
private:
    std::vector<std::vector<Rcpp::IntegerVector> > collected;
};


// [[Rcpp::export(rng=false)]]
Rcpp::IntegerVector fine_tune_label (SEXP Exprs, Rcpp::NumericMatrix scores, Rcpp::List References, 
    double quantile, double tune_thresh, Rcpp::List marker_genes) 
{
    auto mat=beachmat::create_numeric_matrix(Exprs);
    
    matrix_list references; 
    for (size_t i=0; i<References.size(); ++i) {
        references.push_back(beachmat::create_numeric_matrix(References[i]));
    }

    fine_tuner tuner(mat->get_nrow());
    de_markers chooser(marker_genes);

    size_t ncells=mat->get_ncol();
    Rcpp::IntegerVector output(ncells);

    for (size_t c=0; c<ncells; ++c) {
        output[c]=tuner.assign(c, mat.get(), scores, references, quantile, tune_thresh, chooser);
    }

    return output;
}
