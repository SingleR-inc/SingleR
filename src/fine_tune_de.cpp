#include "beachmat/numeric_matrix.h"
#include "fine_tuner.h"
#ifdef _OPENMP
#include <omp.h>
#endif

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
    references.reserve(References.size());
    matrix_list_raw raw_references;
    raw_references.reserve(References.size());

    for (size_t i=0; i<References.size(); ++i) {
        references.push_back(beachmat::create_numeric_matrix(References[i]));
        raw_references.push_back(references.back().get());
    }

    fine_tuner tuner(mat->get_nrow());
    de_markers chooser(marker_genes);

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
            auto tmp=tuner.assign(c, matptr_raw, scores, curref_raw, quantile, tune_thresh, chooser);
            output_id[c]=std::get<0>(tmp);
            output_best[c]=std::get<1>(tmp);
            output_next[c]=std::get<2>(tmp);
        }
    }

    return Rcpp::List::create(output_id, output_best, output_next);
}
