#ifndef SCRAN_PCA_UTILS_GENERAL_HPP
#define SCRAN_PCA_UTILS_GENERAL_HPP

#include <cmath>
#include <algorithm>

namespace scran {

namespace pca_utils {

template<class EigenVector_>
auto process_scale_vector(bool scale, EigenVector_& scale_v) {
    typedef typename EigenVector_::Scalar Scalar;
    if (scale) {
        Scalar total_var = 0;
        for (auto& s : scale_v) {
            if (s) {
                s = std::sqrt(s);
                ++total_var;
            } else {
                s = 1; // avoid division by zero.
            }
        }
        return total_var;
    } else {
        return std::accumulate(scale_v.begin(), scale_v.end(), static_cast<Scalar>(0.0));
    }
}

template<class EigenMatrix_, class EigenVector_>
void clean_up(size_t NC, EigenMatrix_& U, EigenVector_& D) {
    U.array().rowwise() *= D.adjoint().array();
    for (auto& d : D) {
        d = d * d / static_cast<double>(NC - 1);
    }
}

}

}

#endif
