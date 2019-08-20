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
    void operator() (const std::vector<int>& labels, std::set<int>& genes) const {
        genes.clear();
        for (auto l : labels) {
            for (auto l2 : labels) {
                auto& current=collected[l][l2];
                genes.insert(current.begin(), current.end());
            }
        }
        return;
    }
private:
    std::vector<std::vector<Rcpp::IntegerVector> > collected;
};

//' @importFrom Rcpp sourceCpp
//' @useDynLib SingleR
// [[Rcpp::export(rng=false)]]
Rcpp::List fine_tune_label_de (SEXP Exprs, Rcpp::NumericMatrix scores, Rcpp::List References, 
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
    Rcpp::IntegerVector output_id(ncells);
    Rcpp::NumericVector output_best(ncells);
    Rcpp::NumericVector output_next(ncells);

    for (size_t c=0; c<ncells; ++c) {
        auto tmp=tuner.assign(c, mat.get(), scores, references, quantile, tune_thresh, chooser);
        output_id[c]=std::get<0>(tmp);
        output_best[c]=std::get<1>(tmp);
        output_next[c]=std::get<2>(tmp);
    }

    return Rcpp::List::create(output_id, output_best, output_next);
}
