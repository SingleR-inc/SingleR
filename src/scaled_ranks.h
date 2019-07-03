#include <algorithm>
#include <vector>

typedef std::vector<std::pair<double, size_t> > ranked_vector;

template <class IT> 
void scaled_ranks(IT start, IT end, const std::set<int>& chosen, ranked_vector& collected, std::vector<double>& outgoing) {
    size_t slen=chosen.size();

    // Various other bits and pieces.
    collected.reserve(slen);
    collected.clear();
    outgoing.clear();
    outgoing.resize(slen);

    // Sorting all subsetted values (zeroes are handled separately for greater efficiency).
    for (auto i : chosen) {
        const double curval=*(IT + i);
        if (isNA(curval)) { 
            throw std::runtime_error("missing values not supported in quickCluster");
        } else {
            collected.push_back(std::make_pair(curval, s));
        }
    }
    std::sort(collected.begin(), collected.end());

    // Computing tied ranks. 
    double accumulated_rank=0;
    size_t cur_rank=0;
    auto cIt=collected.begin();

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
        o-=center_rank;
        sum_squares+=o*o;
    }

    if (sum_squares==0) {
        for (auto& o : outgoing) {
            o=NA_Real;
        }
    } else {
        sum_squares = std::sqrt(sum_squares)*2;
        for (auto& o : outgoing) {
            o/=sum_squares;
        }
    }
    return;
}
