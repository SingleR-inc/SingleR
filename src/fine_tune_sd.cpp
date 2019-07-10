#include "beachmat/numeric_matrix.h"
#include "fine_tuner.h"

struct sd_markers {
    sd_markers(Rcpp::NumericMatrix Mat, double l) : mat(Mat), limit(l*l) {}
    void operator()(const std::vector<int>& labels, std::set<int>& genes) {
        genes.clear();
        auto it=mat.begin();
        for (int i=0; i<mat.ncol(); ++i, it+=mat.nrow()) {
            // First pass, collect the mean.
            double mean=0;
            for (auto l : labels) {
                mean+=*(it + l);
            }
            mean /= labels.size();

            // Second pass, compute the variance.
            double var=0;
            for (auto l : labels) {
                const double tmp=*(it + l) - mean;
                var+=tmp * tmp;
            }
            var /= labels.size() - 1;

            // Limit is already squared.
            if (var > limit) {
                genes.insert(i);
            }
        }
        return;
    }
private:
    Rcpp::NumericMatrix mat;
    double limit;
};

//' @importFrom Rcpp sourceCpp
//' @useDynLib SingleR
// [[Rcpp::export(rng=false)]]
Rcpp::IntegerVector fine_tune_label_sd (SEXP Exprs, Rcpp::NumericMatrix scores, Rcpp::List References, 
    double quantile, double tune_thresh, Rcpp::NumericMatrix median_mat, double sd_thresh) 
{
    auto mat=beachmat::create_numeric_matrix(Exprs);
    
    matrix_list references; 
    for (size_t i=0; i<References.size(); ++i) {
        references.push_back(beachmat::create_numeric_matrix(References[i]));
    }

    fine_tuner tuner(mat->get_nrow());
    sd_markers chooser(median_mat, sd_thresh);

    size_t ncells=mat->get_ncol();
    Rcpp::IntegerVector output(ncells);

    for (size_t c=0; c<ncells; ++c) {
        output[c]=tuner.assign(c, mat.get(), scores, references, quantile, tune_thresh, chooser);
    }

    return output;
}
