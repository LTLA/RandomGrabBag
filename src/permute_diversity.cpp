#include "Rcpp.h"
#include "Rcpp.h"
#include "pcg_random.hpp"
#include "convert_seed.h"
#include "boost/random.hpp"
#include <algorithm>
#include <vector>
#include <cmath>

template<class ITER>
double compute_discrete_gini(ITER start, ITER end) {
    std::sort(start, end);
    double accumulated=0, current=0;
    for (auto copy=start; copy!=end; ++copy) {
        current += *copy;
        accumulated += current;
    }
    if (current==0) {
        return R_NaReal;
    } else {
        return 1 - accumulated / (current * (current + 1)/2);   
    }
}

template<class ITER>
double compute_hill(ITER start, ITER end, int hill_no) {
    double total=std::accumulate(start, end, 0);
    if (total==0) {
        return R_NaReal;
    }

    double output=0;
    if (hill_no!=1) {
        for (auto copy=start; copy!=end; ++copy) {
            if (*copy) {
                output += std::pow(*copy/total, 1/(1.0-hill_no));
            }
        }
    } else {
        for (auto copy=start; copy!=end; ++copy) {
            if (*copy) {
                const double p=*copy/total;
                output -= p*std::log(p);
            }
            output=std::exp(output);
        }
    }

    return output;
}

//' @useDynLib RepertoireUtils
//' @importFrom Rcpp sourceCpp
//[[Rcpp::export(rng=false)]]
Rcpp::List permute_diversity(Rcpp::IntegerVector left, Rcpp::IntegerVector right, 
    int iterations, bool use_gini, Rcpp::IntegerVector use_hill,
    Rcpp::IntegerVector seed, int stream)
{
    const int ltotal=std::accumulate(left.begin(), left.end(), 0L);
    const int rtotal=std::accumulate(right.begin(), right.end(), 0L);
    const int total=ltotal+rtotal;

    /* Creating the pool by adding the ranked frequencies. This gives a more
     * realistic population composition under the null, in the sense that 
     * the total number of clonotypes is similar to that in the two original
     * groups when the null hypothesis is true. By comparison, literal pooling
     * would double the number of clonotypes and halve the chance of selecting
     * each one, which does not make a lot of sense.
     */
    const auto poolsize=std::max(left.size(), right.size());
    std::vector<int> pool(poolsize);
    {
        std::vector<int> lpool(poolsize);
        std::copy(left.begin(), left.end(), lpool.begin());
        std::copy(right.begin(), right.end(), pool.begin());

        std::sort(lpool.begin(), lpool.end());
        std::sort(pool.begin(), pool.end());

        auto lIt=lpool.begin();
        for (auto pIt=pool.begin(); pIt!=pool.end(); ++pIt, ++lIt) {
            *pIt += *lIt;
        }
    }

    pcg32 rng(dqrng::convert_seed<uint64_t>(seed), stream);
    boost::random::uniform_01<double> runif_gen;

    Rcpp::NumericVector out_gini;
    if (use_gini) { 
        out_gini=Rcpp::NumericVector(iterations);
    }
    Rcpp::NumericMatrix out_hill;
    if (use_hill.size()) {
        out_hill=Rcpp::NumericMatrix(use_hill.size(), iterations);
    }

    std::vector<int> sampled_left, sampled_right;
    sampled_left.reserve(left.size());
    sampled_right.reserve(right.size());

    for (int i=0; i<iterations; ++i) {
        // Sampling to split left/right again.
        int remaining_select=ltotal, remaining_total=total;
        sampled_left.clear();
        sampled_right.clear();

        for (int j=0; j<pool.size(); ++j) {
            int left_assign=0, right_assign=0;
            int curpool=pool[j];

            for (int k=0; k<curpool; ++k) {
                if (remaining_select && static_cast<double>(remaining_select)/remaining_total > runif_gen(rng)) {
                    --remaining_select;
                    ++left_assign;
                } else {
                    ++right_assign;
                }
                --remaining_total;
            }

            if (left_assign) {
                sampled_left.push_back(left_assign);
            }
            if (right_assign) {
                sampled_right.push_back(right_assign);
            }
        }

        // Computing differences in Gini indices and Hill numbers.
        if (use_gini) {
            out_gini[i]=compute_discrete_gini(sampled_left.begin(), sampled_left.end()) - 
                compute_discrete_gini(sampled_right.begin(), sampled_right.end());
        }
        if (use_hill.size()) {
            auto curcol=out_hill.column(i);
            for (auto h=0; h<use_hill.size(); ++h) {
                curcol[h]=compute_hill(sampled_left.begin(), sampled_left.end(), use_hill[h]) -
                    compute_hill(sampled_right.begin(), sampled_right.end(), use_hill[h]);
            }
        }
    }

    return Rcpp::List::create(out_gini, out_hill);
}
