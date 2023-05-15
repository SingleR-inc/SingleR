#ifndef SINGLEPP_CHOOSE_CLASSIC_MARKERS_HPP
#define SINGLEPP_CHOOSE_CLASSIC_MARKERS_HPP

#include "Markers.hpp"
#include "tatami/tatami.hpp"
#include <vector>
#include <cmath>
#include <set>

/**
 * @file ChooseClassicMarkers.hpp
 *
 * @brief Classic method for choosing markers.
 */

namespace singlepp {

/**
 * @brief Classic method for choosing markers.
 *
 * This class implements the classic **SingleR** method for choosing markers from (typically bulk) reference dataasets.
 * We assume that we have a matrix of representative log-expression profiles for each label, typically computed by taking some average across all reference profiles for that label.
 * For the comparison between labels A and B, we define the marker set as the top genes with the largest positive differences in A's profile over B (i.e., the log-fold change).
 * The number of top genes can either be explicitly specified or it can be automatically determined from the number of labels.
 *
 * If multiple references are present, we compute the sum of the log-fold changes across all references with both labels A and B.
 * The ordering of the top genes is then performed using the sum across references.
 * It is assumed that all references have the same number and ordering of features in their rows.
 */
class ChooseClassicMarkers {
public:
    /**
     * @brief Default parameter settings.
     */
    struct Defaults {
        /**
         * See `set_number()` for more details.
         */
        static constexpr int number = -1;

        /**
         * See `set_num_threads()` for more details.
         */
        static constexpr int num_threads = 1;
    };

private:
    int number = Defaults::number;
    int nthreads = Defaults::num_threads;

public:
    /**
     * @param n Number of top genes to use as the marker set in each pairwise comparison.
     * If -1, this is automatically determined from the number of labels, see `number_of_markers()`.
     *
     * @return A reference to this `ChooseClassicMarkers` object.
     */
    ChooseClassicMarkers& set_number(int n = Defaults::number) {
        number = n;
        return *this;
    }

    /**
     * @param nlabels Number of labels in the reference(s).
     *
     * @return An appropriate number of markers for each pairwise comparison.
     *
     * The exact expression is defined as $500 (\frac{2}{3})^{\log_2{N}}$ for $N$ labels,
     * which steadily decreases the markers per comparison as the number of labels increases.
     * This aims to avoid an excessive number of features when dealing with references with many labels.
     */
    static int number_of_markers(int nlabels) {
        return std::round(500.0 * std::pow(2.0/3.0, std::log(static_cast<double>(nlabels)) / std::log(2.0)));
    }

    /**
     * @param n Number of threads to use.
     *
     * @return A reference to this `ChooseClassicMarkers` object.
     */
    ChooseClassicMarkers& set_num_threads(int n = Defaults::num_threads) {
        nthreads = n;
        return *this;
    }

public:
    /**
     * @tparam Matrix A **tatami** matrix.
     * @tparam Label Integer type for the label identity.
     *
     * @param representatives Vector of representative matrices.
     * Each matrix should contain one column per label; each column should have a representative log-expression profile for that label.
     * All matrices should have the same number of rows, corresponding to the same features.
     * @param labels Vector of pointers of length equal to `representatives`.
     * Each array should be of length equal to the number of columns of the corresponding entry of `representatives`.
     * Each value of the array should specify the label for the corresponding column in its matrix.
     * Values should lie in $[0, N)$ for $N$ unique labels across all entries of `labels`.
     *
     * @return A `Markers` object containing the top markers for each pairwise comparison between labels.
     */
    template<class Matrix, typename Label>
    Markers run(const std::vector<const Matrix*>& representatives, const std::vector<const Label*>& labels) const {
        /**
         * @cond
         */
        size_t nrefs = representatives.size();
        if (nrefs != labels.size()) {
            throw std::runtime_error("'representatives' and 'labels' should have the same length");
        }
        if (nrefs == 0) {
            throw std::runtime_error("'representatives' should contain at least one entry");
        }
        size_t ngenes = representatives.front()->nrow();

        // Determining the total number of labels.
        int nlabels = 0;
        for (size_t r = 0; r < nrefs; ++r) {
            size_t ncols = representatives[r]->ncol();
            auto curlab = labels[r];
            for (size_t c = 0; c < ncols; ++c) {
                if (nlabels <= curlab[c]) {
                    nlabels = curlab[c] + 1;
                }
            }
        }

        // Generating mappings.
        std::vector<std::vector<int> > labels_to_index(nrefs, std::vector<int>(nlabels, -1));
        for (size_t r = 0; r < nrefs; ++r) {
            size_t ncols = representatives[r]->ncol();
            auto curlab = labels[r];
            auto& current = labels_to_index[r];
            for (size_t c = 0; c < ncols; ++c) {
                auto& dest = current[curlab[c]];
                if (dest != -1) {
                    throw std::runtime_error("each label should correspond to no more than one column in each reference");
                }
                current[curlab[c]] = c;
            }
        }

        // Generating pairs for compute; this sacrifices some memory for convenience.
        std::vector<std::pair<int, int> > pairs;
        {
            std::set<std::pair<int, int> > pairs0;
            for (size_t r = 0; r < nrefs; ++r) {
                size_t ncols = representatives[r]->ncol();
                auto curlab = labels[r];
                for (size_t c1 = 0; c1 < ncols; ++c1) {
                    for (size_t c2 = 0; c2 < c1; ++c2) {
                        pairs0.emplace(curlab[c1], curlab[c2]);
                    }
                }
            }
            pairs.insert(pairs.end(), pairs0.begin(), pairs0.end());
            std::sort(pairs.begin(), pairs.end());
        }
        size_t npairs = pairs.size();

        Markers output(nlabels, std::vector<std::vector<int> >(nlabels));

        int actual_number = number;
        if (number < 0) {
            actual_number = std::round(500.0 * std::pow(2.0/3.0, std::log(static_cast<double>(nlabels)) / std::log(2.0)));
        } 
        if (actual_number > static_cast<int>(ngenes)) {
            actual_number = ngenes;
        }

#ifndef SINGLEPP_CUSTOM_PARALLEL
        #pragma omp parallel num_threads(nthreads)
        {
#else
        SINGLEPP_CUSTOM_PARALLEL([&](int, size_t start, size_t end) -> void {
#endif
            
            std::vector<std::pair<double, int> > sorter(ngenes), sorted_copy(ngenes);
            typedef typename  Matrix::value_type Value_;
            typedef typename  Matrix::index_type Index_;
            std::vector<Value_> rbuffer(ngenes), lbuffer(ngenes);
            std::vector<std::shared_ptr<tatami::FullDenseExtractor<Value_, Index_> > > rworks(nrefs), lworks(nrefs);

#ifndef SINGLEPP_CUSTOM_PARALLEL
            #pragma omp for
            for (size_t p = 0; p < npairs; ++p) {
#else
            for (size_t p = start; p < end; ++p) {
#endif

                auto curleft = pairs[p].first;
                auto curright = pairs[p].second;

                auto sIt = sorter.begin();
                for (int g = 0; g < ngenes; ++g, ++sIt) {
                    sIt->first = 0;
                    sIt->second = g;
                }

                for (size_t i = 0; i < nrefs; ++i) {
                    const auto& curavail = labels_to_index[i];
                    auto lcol = curavail[curleft];
                    auto rcol = curavail[curright];
                    if (lcol == -1 || rcol == -1) {
                        continue;                            
                    }

                    // Initialize extractors as needed.
                    auto& lwrk = lworks[i];
                    if (!lwrk) {
                        lwrk = representatives[i]->dense_column();
                    }
                    auto lptr = lwrk->fetch(lcol, lbuffer.data());

                    auto& rwrk = rworks[i];
                    if (!rwrk) {
                        rwrk = representatives[i]->dense_column();
                    }
                    auto rptr = rwrk->fetch(rcol, rbuffer.data());

                    auto sIt = sorter.begin();
                    for (int g = 0; g < ngenes; ++g, ++lptr, ++rptr, ++sIt) {
                        sIt->first += *lptr - *rptr; 
                    }
                }

                for (int flip = 0; flip < 2; ++flip) {
                    // Flipping the log-fold changes so we do the reverse. 
                    // We do a copy so that the treatment of tries would be the same as if
                    // we had sorted on the reversed log-fold changes in the first place;
                    // otherwise the first sort might change the order of ties.
                    if (flip) {
                        sorter = sorted_copy;
                        for (auto& s : sorter) {
                            s.first *= -1;
                        }
                    } else {
                        sorted_copy = sorter;
                    }

                    // partial sort is guaranteed to be stable due to the second index resolving ties.
                    std::partial_sort(sorter.begin(), sorter.begin() + actual_number, sorter.end());

                    std::vector<int> stuff;
                    stuff.reserve(actual_number);
                    for (int g = 0; g < actual_number && sorter[g].first < 0; ++g) { // only keeping those with positive log-fold changes (negative after reversing).
                        stuff.push_back(sorter[g].second); 
                    }

                    if (flip) {
                        output[curleft][curright] = std::move(stuff);
                    } else {
                        output[curright][curleft] = std::move(stuff);
                    }
                }
            }

#ifndef SINGLEPP_CUSTOM_PARALLEL
        }
#else    
        }, npairs, nthreads);
#endif        
        /**
         * @endcond
         */

        return output;
    }

    /**
     * @tparam Matrix A **tatami** matrix.
     *
     * @param representative A representative matrix, containing one column per label.
     * Each column should have a representative log-expression profile for that label.
     * @param labels Pointer to an array of length equal to the number of columns in `representative`.
     * Each value of the array should specify the label for the corresponding column. 
     * Values should lie in $[0, N)$ for $N$ unique labels. 
     *
     * @return A `Markers` object containing the top markers for each pairwise comparison between labels.
     */
    template<class Matrix>
    Markers run(const Matrix* representative, const int* labels) const {
        return run(std::vector<const Matrix*>{ representative }, std::vector<const int*>{ labels });
    }
};

}

#endif
