#include "beachmat/numeric_matrix.h"
#include "fine_tuner.h"
#ifdef _OPENMP
#include <omp.h>
#endif

struct sd_markers {
    sd_markers(Rcpp::NumericMatrix Mat, double l) : mat(Mat), limit(l*l) {}
    void operator() (const std::vector<int>& labels, std::set<int>& genes) const {
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
Rcpp::List fine_tune_label_sd (SEXP Exprs, Rcpp::NumericMatrix scores, Rcpp::List References, 
    double quantile, double tune_thresh, Rcpp::NumericMatrix median_mat, double sd_thresh) 
{
    auto mat=beachmat::create_numeric_matrix(Exprs);
    
    matrix_list references;
    references.reserve(References.size());
    matrix_list_raw raw_references;
    raw_references.reserve(References.size());

    for (size_t i=0; i<References.size(); ++i) {
        references.push_back(beachmat::create_numeric_matrix(References[i]));
        raw_references.push_back(references.back().get());
    }

    fine_tuner tuner(mat->get_nrow());
    sd_markers chooser(median_mat, sd_thresh);

    size_t ncells=mat->get_ncol();
    Rcpp::IntegerVector output_id(ncells);
    Rcpp::NumericVector output_best(ncells);
    Rcpp::NumericVector output_next(ncells);

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef _OPENMP
        matrix_list curref;
        matrix_list_raw curref_raw;
        std::unique_ptr<beachmat::numeric_matrix*> matptr=NULL;
        beachmat::numeric_matrix* matptr_raw=NULL;

        #pragma omp critical
        {
            if (omp_get_num_threads()==1 || omp_get_thread_num()==1) {
                curref_raw=raw_references;
                matptr_raw=mat.get();
            } else {
                curref.reserve(references.size());
                curref_raw.reserve(references.size());
                for (size_t i=0; i<curref.size(); ++i) {
                    curref.push_back(references.clone());
                    curref_raw.push_back(curref.back().get());
                }
                matptr=mat.clone();
                matptr_raw=matpr.get();
            }
        }
#else
        matrix_list_raw& curref_raw=raw_references;
        beachmat::numeric_matrix* matptr_raw=mat.get();
#endif

#ifdef _OPENMP
        #pragma omp for static
#endif        
        for (size_t c=0; c<ncells; ++c) {
            auto tmp=tuner.assign(c, matptr_raw, scores, raw_references, quantile, tune_thresh, chooser);
            output_id[c]=std::get<0>(tmp);
            output_best[c]=std::get<1>(tmp);
            output_next[c]=std::get<2>(tmp);
        }
    }

    return Rcpp::List::create(output_id, output_best, output_next);
}
