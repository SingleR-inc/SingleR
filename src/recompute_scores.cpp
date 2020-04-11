#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "compute_scores.h"

#include <vector>
#include <stdexcept>
#include <set>

// [[Rcpp::export(rng=false)]]
Rcpp::RObject recompute_scores(Rcpp::RObject Exprs, Rcpp::IntegerMatrix Labels, 
    Rcpp::List References, Rcpp::List Genes, double quantile) 
{
    auto mat=beachmat::create_numeric_matrix(Exprs);
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();
    Rcpp::NumericVector tmp(ngenes);

    /* Preparing a whole bunch of sanity checks. */

    const size_t nref=References.size();
    matrix_list references; 
    std::vector<int> shift(nref+1);

    for (size_t i=0; i<nref; ++i) {
        Rcpp::List more_references=References[i];
        const size_t nmore=more_references.size();
        shift[i+1]=nmore;
        
        for (size_t j=0; j<nmore; ++j) {
            references.push_back(beachmat::create_numeric_matrix(more_references[j]));
            if (references.back()->get_nrow()!=ngenes) {
                throw std::runtime_error("each entry of 'References' must have number of rows equal to 'Exprs'");
            }
        }
    }

    if (Labels.nrow()!=nref) {
        throw std::runtime_error("'nrow(Labels)' and 'length(References)' must be the same");
    }
    if (Labels.ncol()!=ncells) {
        throw std::runtime_error("'ncol(Labels)' and 'ncol(Exprs)' must be the same");
    }

    if (Genes.size()!=nref) {
        throw std::runtime_error("'Genes' and 'References' must be of the same length");
    }
    std::vector<std::vector<Rcpp::IntegerVector> > genes(nref);
    for (size_t i=0; i<nref; ++i) {
        Rcpp::List current_from=Genes[i];
        const size_t nlab=current_from.size();
        auto& current_to=genes[i];
        current_to.reserve(nlab);

        for (size_t j=0; j<nlab; ++j) {
            current_to.push_back(current_from[j]);
        }
    }

    /* Computing the scores. */

    Rcpp::NumericVector holder_left(ngenes), holder_right(ngenes);
    ranked_vector collected(ngenes);
    std::vector<double> scaled_left(ngenes), scaled_right(ngenes);
    std::vector<double> all_correlations;
    std::vector<int> curlabels(nref);

    Rcpp::NumericMatrix output_best(nref, ncells);

    for (size_t c=0; c<ncells; ++c) {
        mat->get_col(c, holder_left.begin());

        // Finding all of the genes of interest.
        std::set<int> tmp;
        auto all_labels=Labels.column(c);
        for (size_t i=0; i<nref; ++i) {
            Rcpp::IntegerVector current=genes[i][all_labels[i]];
            tmp.insert(current.begin(), current.end());
        }
        std::vector<int> universe(tmp.begin(), tmp.end()); // switch to a more cache-efficient vector.

        // Finding the indices of the labels of interest.
        for (size_t i=0; i<nref; ++i) {
            curlabels[i]=all_labels[i]+shift[i];
        }

        auto current_out=output_best.column(c);
        compute_scores(holder_left, 
            references,
            curlabels, 
            universe,
            quantile,
            holder_right,
            scaled_left,
            scaled_right,
            collected,
            all_correlations,
            current_out.begin()
        );
    }

    return output_best;
}
