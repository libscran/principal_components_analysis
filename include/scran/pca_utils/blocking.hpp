#ifndef SCRAN_PCA_UTILS_MOMENTS_HPP
#define SCRAN_PCA_UTILS_MOMENTS_HPP

#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"
#include "Eigen/Dense"

#include <vector>

#include "general.hpp"
#include "block_weights.hpp"

namespace scran {

namespace pca_utils {

template<typename Index_, class EigenVector_>
struct BlockingDetails {
    std::vector<Index_> block_size;

    bool weighted = false;
    typedef typename EigenVector_::Scalar Weight;

    // The below should only be used if weighted = true.
    std::vector<Weight> per_element_weight;
    Weight total_block_weight = 0;
    EigenVector_ expanded_weights;
};

template<typename Index_, class EigenVector_, typename Block_>
BlockingDetails<Index_, EigenVector_> compute_blocking_details(
    Index_ ncells,
    const Block_* block,
    block_weights::Policy block_weight_policy, 
    const block_weights::VariableParameters& variable_block_weight_parameters) 
{
    auto bsizes = tatami_stats::tabulate_groups(block, ncells);
    BlockingDetails<Index_, EigenVector_> output;
    output.block_size.swap(bsizes);

    if (block_weight_policy == WeightPolicy::NONE) {
        return output;
    }

    output.weighted = true;
    auto& total_weight = output.total_block_weight;
    auto& element_weight = output.per_element_weight;

    for (size_t i = 0; i < nblocks; ++i) {
        auto block_size = output.block_size[i];

        // Computing effective block weights that also incorporate division by the
        // block size. This avoids having to do the division by block size in the
        // 'compute_mean_and_variance_regress()' function.
        if (block_size) {
            typename EigenVector_::Scalar block_weight = 1;
            if (block_weight_policy == WeightPolicy::VARIABLE) {
                block_weight = variable_block_weight(block_size, variable_block_weight_parameters);
            }

            element_weight[i] = block_weight / block_size;
            total_weight += block_weight;
        } else {
            element_weight[i] = 0;
        }
    }

    // Setting a placeholder value to avoid problems with division by zero.
    if (total_weight == 0) {
        total_weight = 1; 
    }

    // Expanding them for multiplication in the IRLBA wrappers.
    auto sqrt_weights = element_weight;
    for (auto& s : sqrt_weights) {
        s = std::sqrt(s);
    }

    auto& expanded = output.expanded_weights;
    expanded.resize(ncells);
    for (size_t i = 0; i < ncells; ++i) {
        expanded.coeffRef(i) = sqrt_weights[block[i]];
    }

    return output;
}

template<typename Num_, typename Value_, typename Index_, typename Block_, typename EigenVector_>
void compute_sparse_mean_and_variance(
    Num_ num_nonzero, 
    const Value_* values, 
    const Index_* indices, 
    const Block* block, 
    const BlockingDetails<Index_, EigenVector_>& block_details,
    Float_* centers,
    Float& variance,
    std::vector<BlockSize_>& block_copy) 
{
    const auto& block_size = block_details.block_size;
    size_t nblocks = block_size.size();

    std::fill_n(center, nblocks, 0);
    for (Num_ i = 0; i < num_nonzero; ++i) {
        centers[block[indices[i]]] += values[i];
    }
    for (size_t b = 0; b < nblocks; ++b) {
        auto bsize = block_size[b];
        if (bsize) {
            centers[b] /= bsize;
        }
    }

    // Computing the variance from the sum of squared differences.
    // This is technically not the correct variance estimate if we
    // were to consider the loss of residual d.f. from estimating
    // the block means, but it's what the PCA sees, so whatever.
    variance = 0;
    std::copy(block_size.begin(), block_size.end(), block_copy.begin());

    if (block_details.weighted) {
        for (Num_ i = 0; i < num_nonzero; ++i) {
            Block_ curb = block[indices[i]];
            auto diff = values[i] - centers[curb];
            variance += diff * diff * block_details.per_element_weight[curb];
            --block_copy[curb];
        }
        for (size_t b = 0; b < nblocks; ++b) {
            auto val = centers[b];
            variance += val * val * block_copy[b] * block_details.per_element_weight[b];
        }
    } else {
        for (Num_ i = 0; i < num_nonzero; ++i) {
            Block_ curb = block[indices[i]];
            auto diff = values[i] - centers[curb];
            variance += diff * diff;
            --block_copy[curb];
        }
        for (size_t b = 0; b < nblocks; ++b) {
            auto val = centers[b];
            variance += val * val * block_copy[b];
        }
    }

    // If we're not dealing with weights, we compute the actual sample
    // variance for easy interpretation (and to match up with the
    // per-PC calculations in pca_utils::clean_up).
    //
    // If we're dealing with weights, the concept of the sample
    // variance becomes somewhat weird, but we just use the same
    // denominator for consistency in pca_utils::clean_up_projected.
    // Magnitude doesn't matter when scaling for
    // pca_utils::process_scale_vector anyway.
    variance /= NR - 1;
}

template<class IrlbaSparseMatrix_, typename Block_, class Index_, class EigenVector_, class EigenMatrix_>
void compute_blockwise_mean_and_variance_realized_sparse(
    const IrlbaSparseMatrix_& emat, // this should be column-major with genes in the columns.
    const Block_* block, 
    const BlockingDetails<Index_, EigenVector_>& block_details,
    EigenMatrix_& centers,
    EigenVector_& variances,
    int nthreads) 
{
    tatami::parallelize([&](size_t, Eigen::Index start, Eigen::Index length) -> void {
        auto NR = emat.rows();
        const auto& values = emat.get_values();
        const auto& indices = emat.get_indices();
        const auto& ptrs = emat.get_pointers();

        size_t nblocks = block_details.block_size.size();
        static_assert(!EigenMatrix_::IsRowMajor);
        auto mptr = centers.data() + static_cast<size_t>(start) * nblocks; // cast to avoid overflow.
        std::vector<Index_> block_copy(nblocks);

        for (size_t c = start, end = start + length; c < end; ++c, mptr += nblocks) {
            auto offset = ptrs[c];
            compute_sparse_mean_and_variance(
                ptrs[c + 1] - offset,
                values.data() + offset,
                indices.data() + offset,
                block,
                block_details,
                mptr,
                variances[c],
                block_copy
            );
        }
    }, emat.cols(), nthreads);
}

template<typename Num_, typename Value_, typename Block_, typename Index_, typename EigenVector_>
void compute_dense_mean_and_variance(
    Num_ number, 
    const Value_* values, 
    const Block* block, 
    const BlockingDetails<Index_, EigenVector_>& block_details,
    Float_* centers,
    Float& variance) 
{
    std::fill_n(centers, nblocks, 0);
    for (size_t r = 0; r < NR; ++r) {
        centers[block[r]] += values[r];
    }
    for (int b = 0; b < nblocks; ++b) {
        const auto& bsize = block_size[b];
        if (bsize) {
            centers[b] /= bsize;
        }
    }

    variance = 0;

    if (block_details.weighted) {
        for (size_t r = 0; r < NR; ++r) {
            auto curb = block[r];
            auto delta = values[r] - centers[curb];
            variance += delta * delta * block_details.per_element_weight[curb];
        }
    } else {
        for (size_t r = 0; r < NR; ++r) {
            auto curb = block[r];
            auto delta = values[r] - centers[curb];
            variance += delta * delta;
        }
    }

    // If we're not dealing with weights, we compute the actual sample
    // variance for easy interpretation (and to match up with the
    // per-PC calculations in pca_utils::clean_up).
    //
    // If we're dealing with weights, the concept of the sample
    // variance becomes somewhat weird, but we just use the same
    // denominator for consistency in pca_utils::clean_up_projected.
    // Magnitude doesn't matter when scaling for
    // pca_utils::process_scale_vector anyway.
    variance/= NR - 1;
}


template<class EigenMatrix_, typename Block_, class Index_, class EigenVector_>
void compute_blockwise_mean_and_variance_realized_dense(
    const EigenMatrix_& emat, // this should be column-major with genes in the columns.
    const Block_* block, 
    const BlockingDetails<Index_, EigenVector_>& block_details,
    EigenMatrix_& centers,
    EigenVector_& variances,
    int nthreads) 
{
    tatami::parallelize([&](size_t, size_t start, size_t length) -> void {
        size_t NR = emat.rows();
        static_assert(!EigenMatrix_::IsRowMajor);
        auto ptr = emat.data() + static_cast<size_t>(start) * NR; 

        const auto& block_size = block_details.block_size;
        size_t nblocks = block_size.size();
        auto mptr = centers.data() + static_cast<size_t>(start) * nblocks; // cast to avoid overflow.

        for (size_t c = start, end = start + length; c < end; ++c, ptr += NR, mptr += nblocks) {
            void compute_dense_mean_and_variance(NR, ptr, block, block_details, mptr, variances[c]);
        }
    }, emat.cols(), nthreads);
}

template<typename Value_, typename Index_, typename Block_, class EigenMatrix_, class EigenVector_>
void compute_blockwise_mean_and_variance_tatami(
    const tatami::Matrix<Value_, Index_>& mat, // this should have genes in the rows!
    const Block_* block, 
    const BlockingDetails<Index_, EigenVector_>& block_details,
    EigenMatrix_& centers,
    EigenVector_& variances,
    int nthreads) 
{
    size_t nblocks = block_details.block_size.size();
    Index_ NR = mat->nrow();
    Index_ NC = mat->ncol();

    if (mat->prefer_rows()) {
        tatami::parallelize([&](size_t, Index_ start, Index_ length) -> void {
            auto mptr = centers.data() + static_cast<size_t>(start) * nblocks; // cast to avoid overflow.
            std::vector<Index_> block_copy(nblocks);

            std::vector<Value_> vbuffer(NC);

            if (mat->is_sparse()) {
                std::vector<Index_> ibuffer(NC);
                auto ext = tatami::consecutive_extractor<true>(mat, true, start, length);
                for (Index_ r = start, end = start + length; r < end; ++r, mptr += nblocks) {
                    auto range = ext->fetch(vbuffer.data(), ibuffer.data());
                    compute_sparse_mean_and_variance(range.number, range.value, range.index, block, block_details, mptr, variances[r], block_copy);
                }
            } else {
                auto ext = tatami::consecutive_extractor<false>(mat, true, start, length);
                for (Index_ r = start, end = start + length; r < end; ++r, mptr += nblocks) {
                    auto ptr = ext->fetch(vbuffer.data());
                    compute_dense_mean_and_variance(NC, ptr, block, block_details, mptr, variances[r]);
                }
            }
        }, NR, nthreads);

    } else {
        typedef typename EigenVector_::Scalar Scalar;

        std::vector<std::pair<size_t, Scalar> > block_multipliers;
        block_multipliers.reserve(nblocks);
        for (size_t b = 0; b < nblocks; ++b) {
            auto bsize = block_size[b];
            if (bsize > 1) { // skipping blocks with NaN variances.
                Scalar mult = bsize - 1; // need to convert variances back into sum of squared differences.
                if (block_details.weighted) {
                    mult *= block_details.per_element_weight[b];
                }
                multipliers.emplace_back(b, mult);
            }
        }

        tatami::parallelize([&](size_t, Index_ start, Index_ length) -> void {
            std::vector<std::vector<Scalar_> > re_centers, re_variances;
            re_centers.reserve(nblocks);
            re_variances.reserve(nblocks);
            for (size_t b = 0; b < nblocks; ++b) {
                re_centers.emplace_back(length);
                re_variances.emplace_back(length);
            }

            std::vector<Value_> vbuffer(length);

            if (mat->is_sparse()) {
                std::vector<tatami_stats::variances::RunningSparse<Scalar_, Value_, Index_> > running;
                running.reserve(nblocks);
                for (size_t b = 0; b < nblocks; ++b) {
                    running.emplace_back(re_centers[b].data(), re_variances[b].data(), /* skip_nan = */ false, /* subtract = */ start);
                }

                std::vector<Index_> ibuffer(length);
                auto ext = tatami::consecutive_extractor<true>(mat, 0, NC, false, start, length);
                for (Index_ c = 0; c < NC; ++c) {
                    auto range = ext->fetch(vbuffer.data(), ibuffer.data());
                    running[block[c]].add(range.value, range.index, range.number);
                }

                for (size_t b = 0; b < nblocks; ++b) {
                    running[b].finish();
                }

            } else {
                std::vector<tatami_stats::variances::RunningDense<Scalar_, Value_, Index_> > running;
                running.reserve(nblocks);
                for (size_t b = 0; b < nblocks; ++b) {
                    running.emplace_back(re_centers[b].data(), re_variances[b].data(), /* skip_nan = */ false);
                }

                auto ext = tatami::consecutive_extractor<false>(mat, 0, NC, false, start, length);
                for (Index_ c = 0; c < NC; ++c) {
                    auto ptr = ext->fetch(vbuffer.data());
                    running[block[c]].add(ptr);
                }

                for (size_t b = 0; b < nblocks; ++b) {
                    running[b].finish();
                }
            }

            auto mptr = centers.data() + static_cast<size_t>(start) * nblocks; // cast to avoid overflow.
            for (Index_ r = 0; r < length; ++r, mptr += nblocks) {
                for (size_t b = 0; b < nblocks; ++b) {
                    mptr[b] = re_centers[b][r];
                }
            }

            for (Index_ r = 0; r < length; ++r, mptr += nblocks) {
                for (const auto& bm : block_multipliers) {
                    my_var += re_variances[bm.first][r] * bm.second;
                }

                // If we're not dealing with weights, we compute the actual sample
                // variance for easy interpretation (and to match up with the
                // per-PC calculations in pca_utils::clean_up).
                //
                // If we're dealing with weights, the concept of the sample
                // variance becomes somewhat weird, but we just use the same
                // denominator for consistency in pca_utils::clean_up_projected.
                // Magnitude doesn't matter when scaling for
                // pca_utils::process_scale_vector anyway.
                my_var /= NC - 1;
            }
        }, NR, nthreads);
    }
}

inline void project_sparse_matrix(const SparseMatrix& emat, Eigen::MatrixXd& pcs, const Eigen::MatrixXd& rotation, bool scale, const Eigen::VectorXd& scale_v, int nthreads) {
    size_t nvec = rotation.cols();
    size_t nrow = emat.rows();
    size_t ncol = emat.cols();

    pcs.resize(nvec, nrow); // used a transposed version for more cache efficiency.
    pcs.setZero();

    const auto& x = emat.get_values();
    const auto& i = emat.get_indices();
    const auto& p = emat.get_pointers();

    if (nthreads == 1) {
        Eigen::VectorXd multipliers(nvec);
        for (size_t c = 0; c < ncol; ++c) {
            multipliers.noalias() = rotation.row(c);
            if (scale) {
                multipliers.array() /= scale_v[c];
            }

            auto start = p[c], end = p[c + 1];
            for (size_t s = start; s < end; ++s) {
                pcs.col(i[s]).noalias() += x[s] * multipliers;
            }
        }
    } else {
        const auto& row_nonzero_starts = emat.get_secondary_nonzero_starts();

#ifndef IRLBA_CUSTOM_PARALLEL
        #pragma omp parallel for num_threads(nthreads)
        for (size_t t = 0; t < nthreads; ++t) {
#else
        IRLBA_CUSTOM_PARALLEL(nthreads, [&](size_t t) -> void { 
#endif

            const auto& starts = row_nonzero_starts[t];
            const auto& ends = row_nonzero_starts[t + 1];
            Eigen::VectorXd multipliers(nvec);

            for (size_t c = 0; c < ncol; ++c) {
                multipliers.noalias() = rotation.row(c);
                if (scale) {
                    multipliers.array() /= scale_v[c];
                }

                auto start = starts[c], end = ends[c];
                for (size_t s = start; s < end; ++s) {
                    pcs.col(i[s]).noalias() += x[s] * multipliers;
                }
            }

#ifndef IRLBA_CUSTOM_PARALLEL
        }
#else
        });
#endif
    }
}

template<bool rows_are_dims_>
void clean_up_projected(Eigen::MatrixXd& proj, Eigen::VectorXd& D) {
    // Empirically centering to give nice centered PCs, because we can't
    // guarantee that the projection is centered in this manner.
    if constexpr(rows_are_dims_) {
        for (size_t i = 0, iend = proj.rows(); i < iend; ++i) {
            proj.row(i).array() -= proj.row(i).sum() / proj.cols();
        }
    } else {
        for (size_t i = 0, iend = proj.cols(); i < iend; ++i) {
            proj.col(i).array() -= proj.col(i).sum() / proj.rows();
        }
    }

    // Just dividing by the number of observations - 1 regardless of weighting.
    double denom = (rows_are_dims_ ? proj.cols() : proj.rows()) - 1;
    for (auto& d : D) {
        d = d * d / denom;
    }
}

}

}

#endif
