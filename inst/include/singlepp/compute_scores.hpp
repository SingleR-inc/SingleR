#ifndef SINGLEPP_COMPUTE_SCORES_HPP
#define SINGLEPP_COMPUTE_SCORES_HPP

#include "macros.hpp"

#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

namespace singlepp {

inline double correlations_to_scores (std::vector<double>& correlations, double quantile) {
    const size_t ncells=correlations.size();
    if (ncells==0) {
        return std::numeric_limits<double>::quiet_NaN();
    } else if (quantile==1 || ncells==1) {
        return *std::max_element(correlations.begin(), correlations.end());
    } else {
        const double denom = ncells - 1; 
        const double prod = denom * quantile;
        const size_t left = std::floor(prod);
        const size_t right = std::ceil(prod);

        std::nth_element(correlations.begin(), correlations.begin() + right, correlations.end());
        const double rightval=correlations[right];
        if (right == left) {
            return rightval;
        }

        // Do NOT be tempted to do the second nth_element with the end at begin()+right;
        // this does not handle ties properly.
        std::nth_element(correlations.begin(), correlations.begin() + left, correlations.end());
        const double leftval=correlations[left];

        // `quantile - left / denom` represents the gap to the smaller quantile,
        // while `right / denom - quantile` represents the gap from the larger quantile.
        // The size of the gap is used as the weight for the _other_ quantile, i.e., 
        // the closer you are to a quantile, the higher the weight.
        // We convert these into proportions by dividing by their sum, i.e., `1/denom`.
        const double leftweight = right - prod;
        const double rightweight = prod - left;

        return rightval * rightweight + leftval * leftweight;
    }
}

template<class Vector>
double distance_to_correlation(size_t n, const Vector& p1, const Vector& p2) {
    double d2 = 0;
    for (size_t i = 0; i < n; ++i) {
        double tmp = p1[i] - p2[i];
        d2 += tmp * tmp;
    }
    return 1 - 2 * d2;
}

}

#endif
