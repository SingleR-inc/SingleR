#include "Rcpp.h"
#include <algorithm>
#include <vector>

typedef std::vector<std::pair<double, size_t> > ranked_vector;

template <class IT> 
void scaled_ranks(IT start, const std::vector<int>& chosen, ranked_vector& collected, std::vector<double>& outgoing, bool na_aware = false) {
    size_t slen=chosen.size();
    collected.clear();
    collected.reserve(slen);

    // Sorting all subsetted values (zeroes are handled separately for greater efficiency).
    size_t s=0;
    for (auto i : chosen) {
        const double curval=*(start + i);
        if (ISNA(curval)) { 
            if (!na_aware) {
                throw std::runtime_error("missing values not supported in SingleR");
            }
        } else {
            collected.push_back(std::make_pair(curval, s));
        }
        ++s;
    }
    std::sort(collected.begin(), collected.end());

    // Computing tied ranks. 
    size_t cur_rank=0;
    auto cIt=collected.begin();
    outgoing.resize(slen, na_aware ? R_NaReal : 0);

    while (cIt!=collected.end()) {
        auto copy=cIt;
        ++copy;
        double accumulated_rank=cur_rank;
        ++cur_rank;

        while (copy!=collected.end() && copy->first==cIt->first) {
            accumulated_rank+=cur_rank;
            ++cur_rank;
            ++copy;
        }

        double mean_rank=accumulated_rank/(copy-cIt);
        while (cIt!=copy) {
            outgoing[cIt->second]=mean_rank;
            ++cIt;
        }
    }

    // Mean-adjusting and converting to cosine values.
    double sum_squares=0;
    const double center_rank=static_cast<double>(slen-1)/2;
    for (auto& o : outgoing) {
        if (na_aware && ISNA(o)) {
            continue;
        }
        o-=center_rank;
        sum_squares+=o*o;
    }

    // Special behaviour for no-variance genes; these are left as all-zero scaled ranks.
    sum_squares = std::max(sum_squares, 0.00000001);
    sum_squares = std::sqrt(sum_squares)*2;
    for (auto& o : outgoing) {
        if (na_aware && ISNA(o)) {
            continue;
        }
        o/=sum_squares;
    }

    return;
}
