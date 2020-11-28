#include "Rcpp.h"
#include "beachmat3/beachmat.h"
#include "compute_scores.h"

#include <vector>
#include <stdexcept>
#include <set>

void preflight(std::vector<matrix_list>& references, 
    std::vector<std::vector<Rcpp::IntegerVector> >& genes,
    Rcpp::List& References, Rcpp::IntegerMatrix& Labels, Rcpp::List& Genes,
    size_t ngenes, size_t ncells)  
{
    const size_t nref=References.size();
    for (size_t i=0; i<nref; ++i) {
        Rcpp::List more_references=References[i];
        const size_t nmore=more_references.size();
        auto& currefs=references[i];

        for (size_t j=0; j<nmore; ++j) {
            currefs.push_back(beachmat::read_lin_block(more_references[j]));
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
    for (size_t i=0; i<nref; ++i) {
        Rcpp::List current_from=Genes[i];
        const size_t nlab=current_from.size();
        auto& current_to=genes[i];
        current_to.reserve(nlab);

        for (size_t j=0; j<nlab; ++j) {
            current_to.push_back(current_from[j]);
        }
    }

    return;
}

std::vector<int> identify_genes (Rcpp::IntegerMatrix::Column& all_labels, 
    const std::vector<std::vector<Rcpp::IntegerVector> >& genes,
    size_t nref)
{
    std::set<int> tmp;
    for (size_t i=0; i<nref; ++i) {
        Rcpp::IntegerVector current=genes[i][all_labels[i]];
        tmp.insert(current.begin(), current.end());
    }
    return std::vector<int>(tmp.begin(), tmp.end()); // switch to a more cache-efficient vector.
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject recompute_scores(Rcpp::List Groups, 
    Rcpp::RObject Exprs, Rcpp::IntegerMatrix Labels, 
    Rcpp::List References, Rcpp::List Genes, double quantile) 
{
    auto mat = beachmat::read_lin_block(Exprs);
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();

    /* Preparing a whole bunch of sanity checks. */
    const size_t nref=References.size();
    const size_t ngroups=Groups.size();
    std::vector<matrix_list> references(nref);
    std::vector<std::vector<Rcpp::IntegerVector> > genes(nref);
    preflight(references, genes, References, Labels, Genes, ngenes, ncells);

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
        auto all_labels=Labels.column(curgroup[0]);
        auto universe=identify_genes(all_labels, genes, nref);

        // Computing the references.
        for (size_t i=0; i<nref; ++i) {
            auto current=references[i][all_labels[i]].get();
            const size_t curncells=current->get_ncol();

            auto& scaled_right_set=scaled_rights[i];
            scaled_right_set.resize(curncells);

            for (size_t c=0; c<curncells; ++c) {
                auto ptr = current->get_col(c, holder_right.begin());
                scaled_ranks(ptr, universe, collected, scaled_right_set[c]);
            }
        }

        // Looping through the cells and computing the scores. 
        for (size_t t=0; t<cursize; ++t) {
            auto ptr = mat->get_col(curgroup[t], holder_left.begin());
            scaled_ranks(ptr, universe, collected, scaled_left);
            auto outcol=output.column(curgroup[t]);

            for (size_t i=0; i<nref; ++i) {
                const auto& cur_set=scaled_rights[i];
                all_correlations.clear();
                for (size_t c=0; c<cur_set.size(); ++c) {
                    const auto& curright=cur_set[c];

                    double dist = 0;
                    auto crIt = curright.begin();
                    for (auto slIt = scaled_left.begin(); slIt != scaled_left.end(); ++slIt, ++crIt) {
                        const double tmp = (*slIt) - (*crIt);
                        dist += tmp*tmp;
                    }
                    all_correlations.push_back(1 - 2*dist);
                }

                outcol[i]=correlations_to_scores(all_correlations, quantile);
            }
        }
    }

    return output;
}

/* A slightly less efficient version that is NA aware. */

// [[Rcpp::export(rng=false)]]
Rcpp::RObject recompute_scores_with_na(Rcpp::List Groups, 
    Rcpp::RObject Exprs, Rcpp::IntegerMatrix Labels, 
    Rcpp::List References, Rcpp::List Genes, double quantile) 
{
    auto mat = beachmat::read_lin_block(Exprs);
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();

    /* Preparing a whole bunch of sanity checks. */
    const size_t nref=References.size();
    const size_t ngroups=Groups.size();
    std::vector<matrix_list> references(nref);
    std::vector<std::vector<Rcpp::IntegerVector> > genes(nref);
    preflight(references, genes, References, Labels, Genes, ngenes, ncells);

    /* Computing the scores. */
    Rcpp::NumericVector holder_left(ngenes), holder_right(ngenes);
    std::vector<double> copy_left(ngenes);
    ranked_vector collected(ngenes);
    std::vector<double> scaled_left(ngenes), scaled_right(ngenes);
    std::vector<double> all_correlations;

    Rcpp::NumericMatrix output(nref, ncells);

    for (size_t g=0; g<ngroups; ++g) {
        Rcpp::IntegerVector curgroup=Groups[g];
        const size_t cursize=curgroup.size();
        if (cursize==0) { 
            continue;
        }

        // Finding all of the genes of interest.
        auto all_labels=Labels.column(curgroup[0]);
        auto universe=identify_genes(all_labels, genes, nref);

        // Looping through the cells and computing the scores. 
        for (size_t t=0; t<cursize; ++t) {
            auto lptr = mat->get_col(curgroup[t], holder_left.begin());
            auto outcol=output.column(curgroup[t]);

            for (size_t i=0; i<nref; ++i) {
                auto current=references[i][all_labels[i]].get();
                const size_t curncells=current->get_ncol();
                std::copy(lptr, lptr + ngenes, copy_left.begin());
                all_correlations.clear();

                for (size_t c=0; c<curncells; ++c) {
                    auto rptr = current->get_col(c, holder_right.begin());
                    scaled_ranks(rptr, universe, collected, scaled_right, true);

                    if (c==0) {
                        // Infect copy_left with the same pattern of NA's,
                        // under the assumption that the reference's NAs are
                        // the same for all its cells.
                        for (auto& u : universe) {
                            if (ISNA(rptr[u])) {
                                copy_left[u]=R_NaReal;
                            }
                        }
                        scaled_ranks(copy_left.begin(), universe, collected, scaled_left, true);
                    }

                    double dist=0;
                    auto srIt=scaled_right.begin();
                    for (auto slIt=scaled_left.begin(); slIt != scaled_left.end(); ++slIt, ++srIt) {
                        if (!ISNA(*slIt)) {
                            const double tmp = (*slIt) - (*srIt);
                            dist+=tmp*tmp;
                        }
                    }
                    all_correlations.push_back(1 - 2*dist);
                }

                outcol[i]=correlations_to_scores(all_correlations, quantile);
            }
        }
    }

    return output;
}
