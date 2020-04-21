#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "compute_scores.h"

#include <vector>
#include <stdexcept>
#include <set>

// [[Rcpp::export(rng=false)]]
Rcpp::RObject recompute_scores(Rcpp::List Groups, 
    Rcpp::RObject Exprs, Rcpp::IntegerMatrix Labels, 
    Rcpp::List References, Rcpp::List Genes, double quantile) 
{
    auto mat=beachmat::create_numeric_matrix(Exprs);
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();
    Rcpp::NumericVector tmp(ngenes);

    /* Preparing a whole bunch of sanity checks. */

    const size_t ngroups=Groups.size();
    const size_t nref=References.size();
    std::vector<matrix_list> references(nref);

    for (size_t i=0; i<nref; ++i) {
        Rcpp::List more_references=References[i];
        const size_t nmore=more_references.size();
        auto& currefs=references[i];

        for (size_t j=0; j<nmore; ++j) {
            currefs.push_back(beachmat::create_numeric_matrix(more_references[j]));
            if (currefs.back()->get_nrow()!=ngenes) {
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
    std::vector<double> scaled_left(ngenes);
    std::vector<std::vector<std::vector<double> > > scaled_rights(nref);
    std::vector<double> all_correlations;

    Rcpp::NumericMatrix output(nref, ncells);

    for (size_t g=0; g<ngroups; ++g) {
        Rcpp::IntegerVector curgroup=Groups[g];
        const size_t cursize=curgroup.size();
        if (cursize==0) { 
            continue;
        }

        // Finding all of the genes of interest.
        std::set<int> tmp;
        auto all_labels=Labels.column(curgroup[0]);
        for (size_t i=0; i<nref; ++i) {
            Rcpp::IntegerVector current=genes[i][all_labels[i]];
            tmp.insert(current.begin(), current.end());
        }
        std::vector<int> universe(tmp.begin(), tmp.end()); // switch to a more cache-efficient vector.

        // Computing the references.
        for (size_t i=0; i<nref; ++i) {
            auto current=references[i][all_labels[i]].get();
            const size_t ncells=current->get_ncol();

            auto& scaled_right_set=scaled_rights[i];
            scaled_right_set.resize(ncells);

            for (size_t c=0; c<ncells; ++c) {
                current->get_col(c, holder_right.begin());
                scaled_ranks(holder_right.begin(), universe, collected, scaled_right_set[c]);
            }
        }

        // Looping through the cells and computing the scores. 
        for (size_t t=0; t<cursize; ++t) {
            mat->get_col(curgroup[t], holder_left.begin());
            scaled_ranks(holder_left.begin(), universe, collected, scaled_left);
            auto outcol=output.column(curgroup[t]);

            for (size_t i=0; i<nref; ++i) {
                const auto& cur_set=scaled_rights[i];
                all_correlations.clear();
                for (size_t c=0; c<cur_set.size(); ++c) {
                    const auto& curright=cur_set[c];

                    double dist=0;
                    for (size_t j=0; j<scaled_left.size(); ++j) {
                        const double tmp=scaled_left[j] - curright[j];
                        dist+=tmp*tmp;
                    }
                    all_correlations.push_back(1 - 2*dist);
                }

                outcol[i]=correlations_to_scores(all_correlations, quantile);
            }
        }
    }

    return output;
}
