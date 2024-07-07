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

template<class EigenVector_, typename Index_, typename Block_>
BlockingDetails<Index_, EigenVector_> compute_blocking_details(
    Index_ ncells,
    const Block_* block,
    block_weights::Policy block_weight_policy, 
    const block_weights::VariableParameters& variable_block_weight_parameters) 
{
    auto bsizes = tatami_stats::tabulate_groups(block, ncells);
    BlockingDetails<Index_, EigenVector_> output;
    output.block_size.swap(bsizes);

    if (block_weight_policy == block_weights::Policy::NONE) {
        return output;
    }

    size_t nblocks = bsizes.size();
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
            if (block_weight_policy == block_weights::Policy::VARIABLE) {
                block_weight = block_weights::compute_variable(block_size, variable_block_weight_parameters);
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
    for (Index_ i = 0; i < ncells; ++i) {
        expanded.coeffRef(i) = sqrt_weights[block[i]];
    }

    return output;
}

/*****************************************************************
 ************ Computing the blockwise mean and variance **********
 *****************************************************************/

template<typename Num_, typename Value_, typename Index_, typename Block_, typename EigenVector_, typename Float_>
void compute_sparse_mean_and_variance(
    Num_ num_nonzero, 
    const Value_* values, 
    const Index_* indices, 
    const Block_* block, 
    const BlockingDetails<Index_, EigenVector_>& block_details,
    Float_* centers,
    Float_& variance,
    std::vector<Index_>& block_copy,
    Num_ num_all)
{
    const auto& block_size = block_details.block_size;
    size_t nblocks = block_size.size();

    std::fill_n(centers, nblocks, 0);
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

    // COMMENT ON DENOMINATOR:
    // If we're not dealing with weights, we compute the actual sample
    // variance for easy interpretation (and to match up with the
    // per-PC calculations in pca_utils::clean_up).
    //
    // If we're dealing with weights, the concept of the sample
    // variance becomes somewhat weird, but we just use the same
    // denominator for consistency in pca_utils::clean_up_projected.
    // Magnitude doesn't matter when scaling for
    // pca_utils::process_scale_vector anyway.
    variance /= num_all - 1;
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
        size_t ncells = emat.rows();
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
                block_copy,
                ncells
            );
        }
    }, emat.cols(), nthreads);
}

template<typename Num_, typename Value_, typename Block_, typename Index_, typename EigenVector_, typename Float_>
void compute_dense_mean_and_variance(
    Num_ number, 
    const Value_* values, 
    const Block_* block, 
    const BlockingDetails<Index_, EigenVector_>& block_details,
    Float_* centers,
    Float_& variance) 
{
    const auto& block_size = block_details.block_size;
    size_t nblocks = block_size.size();
    std::fill_n(centers, nblocks, 0);
    for (Num_ r = 0; r < number; ++r) {
        centers[block[r]] += values[r];
    }
    for (size_t b = 0; b < nblocks; ++b) {
        const auto& bsize = block_size[b];
        if (bsize) {
            centers[b] /= bsize;
        }
    }

    variance = 0;

    if (block_details.weighted) {
        for (Num_ r = 0; r < number; ++r) {
            auto curb = block[r];
            auto delta = values[r] - centers[curb];
            variance += delta * delta * block_details.per_element_weight[curb];
        }
    } else {
        for (Num_ r = 0; r < number; ++r) {
            auto curb = block[r];
            auto delta = values[r] - centers[curb];
            variance += delta * delta;
        }
    }

    variance /= number - 1; // See COMMENT ON DENOMINATOR above.
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
        size_t ncells = emat.rows();
        static_assert(!EigenMatrix_::IsRowMajor);
        auto ptr = emat.data() + static_cast<size_t>(start) * ncells; 

        const auto& block_size = block_details.block_size;
        size_t nblocks = block_size.size();
        auto mptr = centers.data() + static_cast<size_t>(start) * nblocks; // cast to avoid overflow.

        for (size_t c = start, end = start + length; c < end; ++c, ptr += ncells, mptr += nblocks) {
            compute_dense_mean_and_variance(ncells, ptr, block, block_details, mptr, variances[c]);
        }
    }, emat.cols(), nthreads);
}

template<typename Value_, typename Index_, typename Block_, class EigenMatrix_, class EigenVector_>
void compute_blockwise_mean_and_variance_tatami(
    const tatami::Matrix<Value_, Index_>* mat, // this should have genes in the rows!
    const Block_* block, 
    const BlockingDetails<Index_, EigenVector_>& block_details,
    EigenMatrix_& centers,
    EigenVector_& variances,
    int nthreads) 
{
    const auto& block_size = block_details.block_size;
    size_t nblocks = block_size.size();
    Index_ ngenes = mat->nrow();
    Index_ ncells = mat->ncol();

    if (mat->prefer_rows()) {
        tatami::parallelize([&](size_t, Index_ start, Index_ length) -> void {
            static_assert(!EigenMatrix_::IsRowMajor);
            auto mptr = centers.data() + static_cast<size_t>(start) * nblocks; // cast to avoid overflow.
            std::vector<Index_> block_copy(nblocks);

            std::vector<Value_> vbuffer(ncells);

            if (mat->is_sparse()) {
                std::vector<Index_> ibuffer(ncells);
                auto ext = tatami::consecutive_extractor<true>(mat, true, start, length);
                for (Index_ r = start, end = start + length; r < end; ++r, mptr += nblocks) {
                    auto range = ext->fetch(vbuffer.data(), ibuffer.data());
                    compute_sparse_mean_and_variance(range.number, range.value, range.index, block, block_details, mptr, variances[r], block_copy, ncells);
                }
            } else {
                auto ext = tatami::consecutive_extractor<false>(mat, true, start, length);
                for (Index_ r = start, end = start + length; r < end; ++r, mptr += nblocks) {
                    auto ptr = ext->fetch(vbuffer.data());
                    compute_dense_mean_and_variance(ncells, ptr, block, block_details, mptr, variances[r]);
                }
            }
        }, ngenes, nthreads);

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
                block_multipliers.emplace_back(b, mult);
            }
        }

        tatami::parallelize([&](size_t, Index_ start, Index_ length) -> void {
            std::vector<std::vector<Scalar> > re_centers, re_variances;
            re_centers.reserve(nblocks);
            re_variances.reserve(nblocks);
            for (size_t b = 0; b < nblocks; ++b) {
                re_centers.emplace_back(length);
                re_variances.emplace_back(length);
            }

            std::vector<Value_> vbuffer(length);

            if (mat->is_sparse()) {
                std::vector<tatami_stats::variances::RunningSparse<Scalar, Value_, Index_> > running;
                running.reserve(nblocks);
                for (size_t b = 0; b < nblocks; ++b) {
                    running.emplace_back(length, re_centers[b].data(), re_variances[b].data(), /* skip_nan = */ false, /* subtract = */ start);
                }

                std::vector<Index_> ibuffer(length);
                auto ext = tatami::consecutive_extractor<true>(mat, false, static_cast<Index_>(0), ncells, start, length);
                for (Index_ c = 0; c < ncells; ++c) {
                    auto range = ext->fetch(vbuffer.data(), ibuffer.data());
                    running[block[c]].add(range.value, range.index, range.number);
                }

                for (size_t b = 0; b < nblocks; ++b) {
                    running[b].finish();
                }

            } else {
                std::vector<tatami_stats::variances::RunningDense<Scalar, Value_, Index_> > running;
                running.reserve(nblocks);
                for (size_t b = 0; b < nblocks; ++b) {
                    running.emplace_back(length, re_centers[b].data(), re_variances[b].data(), /* skip_nan = */ false);
                }

                auto ext = tatami::consecutive_extractor<false>(mat, false, static_cast<Index_>(0), ncells, start, length);
                for (Index_ c = 0; c < ncells; ++c) {
                    auto ptr = ext->fetch(vbuffer.data());
                    running[block[c]].add(ptr);
                }

                for (size_t b = 0; b < nblocks; ++b) {
                    running[b].finish();
                }
            }

            static_assert(!EigenMatrix_::IsRowMajor);
            auto mptr = centers.data() + static_cast<size_t>(start) * nblocks; // cast to avoid overflow.
            for (Index_ r = 0; r < length; ++r, mptr += nblocks) {
                for (size_t b = 0; b < nblocks; ++b) {
                    mptr[b] = re_centers[b][r];
                }

                auto& my_var = variances[start + r];
                my_var = 0;
                for (const auto& bm : block_multipliers) {
                    my_var += re_variances[bm.first][r] * bm.second;
                }
                my_var /= ncells - 1; // See COMMENT ON DENOMINATOR above.
            }
        }, ngenes, nthreads);
    }
}

/******************************************************************
 ************ Project matrices on their rotation vectors **********
 ******************************************************************/

template<class EigenMatrix_, class EigenVector_>
const EigenMatrix_& scale_rotation_matrix(const EigenMatrix_& rotation, bool scale, const EigenVector_& scale_v, EigenMatrix_& tmp) {
    if (scale) {
        tmp = (rotation.array().colwise() / scale_v.array()).matrix();
        return tmp;
    } else {
        return rotation;
    }
}

template<class IrlbaSparseMatrix_, class EigenMatrix_>
inline void project_matrix_realized_sparse(
    const IrlbaSparseMatrix_& emat, // cell in rows, genes in the columns, CSC.
    EigenMatrix_& components, // dims in rows, cells in columns
    const EigenMatrix_& scaled_rotation, // genes in rows, dims in columns
    int nthreads) 
{
    size_t rank = scaled_rotation.cols();
    size_t ncells = emat.rows();
    size_t ngenes = emat.cols();

    components.resize(rank, ncells); // used a transposed version for more cache efficiency.
    components.setZero();

    const auto& x = emat.get_values();
    const auto& i = emat.get_indices();
    const auto& p = emat.get_pointers();

    if (nthreads == 1) {
        Eigen::VectorXd multipliers(rank);
        for (size_t g = 0; g < ngenes; ++g) {
            multipliers.noalias() = scaled_rotation.row(g);
            auto start = p[g], end = p[g + 1];
            for (size_t s = start; s < end; ++s) {
                components.col(i[s]).noalias() += x[s] * multipliers;
            }
        }

    } else {
        const auto& row_nonzero_starts = emat.get_secondary_nonzero_starts();

#ifndef IRLBA_CUSTOM_PARALLEL
#ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads)
#endif
        for (int t = 0; t < nthreads; ++t) {
#else
        IRLBA_CUSTOM_PARALLEL(nthreads, [&](size_t t) -> void { 
#endif

            const auto& starts = row_nonzero_starts[t];
            const auto& ends = row_nonzero_starts[t + 1];
            Eigen::VectorXd multipliers(rank);

            for (size_t g = 0; g < ngenes; ++g) {
                multipliers.noalias() = scaled_rotation.row(g);
                auto start = starts[g], end = ends[g];
                for (size_t s = start; s < end; ++s) {
                    components.col(i[s]).noalias() += x[s] * multipliers;
                }
            }

#ifndef IRLBA_CUSTOM_PARALLEL
        }
#else
        });
#endif
    }
}

template<typename Value_, typename Index_, class EigenMatrix_>
inline void project_matrix_transposed_tatami(
    const tatami::Matrix<Value_, Index_>* mat, // genes in rows, cells in columns
    EigenMatrix_& components, // dims in rows, cells in columns
    const EigenMatrix_& scaled_rotation, // genes in rows, dims in columns
    int nthreads) 
{
    size_t rank = scaled_rotation.cols();
    auto ngenes = mat->nrow(), ncells = mat->ncol();
    typedef typename EigenMatrix_::Scalar Scalar;

    components.resize(rank, ncells); // used a transposed version for more cache efficiency.

    if (mat->prefer_rows()) {
        tatami::parallelize([&](size_t, Index_ start, Index_ length) -> void {
            static_assert(!EigenMatrix_::IsRowMajor);
            auto vptr = scaled_rotation.data();
            std::vector<Value_> vbuffer(length);

            std::vector<std::vector<Scalar> > local_buffers; // create separate buffers to avoid false sharing.
            local_buffers.reserve(rank);
            for (size_t r = 0; r < rank; ++r) {
                local_buffers.emplace_back(length);
            }

            if (mat->is_sparse()) {
                std::vector<Index_> ibuffer(length);
                auto ext = tatami::consecutive_extractor<true>(mat, true, static_cast<Index_>(0), ngenes, start, length);
                for (Index_ g = 0; g < ngenes; ++g) {
                    auto range = ext->fetch(vbuffer.data(), ibuffer.data());
                    for (size_t r = 0; r < rank; ++r) {
                        auto mult = vptr[r + static_cast<size_t>(ngenes) + static_cast<size_t>(g)]; // cast to avoid overflow.
                        auto& local_buffer = local_buffers[r];
                        for (Index_ i = 0; i < range.number; ++i) {
                            auto c = range.index[i];
                            local_buffer[c] += range.value[i] * mult;
                        }
                    }
                }

            } else {
                auto ext = tatami::consecutive_extractor<false>(mat, true, static_cast<Index_>(0), ngenes, start, length);
                for (Index_ g = 0; g < ngenes; ++g) {
                    auto ptr = ext->fetch(vbuffer.data());
                    for (size_t r = 0; r < rank; ++r) {
                        auto mult = vptr[r + static_cast<size_t>(ngenes) + static_cast<size_t>(g)]; // cast to avoid overflow.
                        auto& local_buffer = local_buffers[r];
                        for (Index_ c = 0; c < length; ++c) {
                            local_buffer[c] += ptr[c] * mult;
                        }
                    }
                }
            }

            for (size_t r = 0; r < rank; ++r) {
                for (Index_ c = 0; c < length; ++c) {
                    components.coeffRef(r, c + start) = local_buffers[r][c];
                }
            }

        }, ncells, nthreads);

    } else {
        tatami::parallelize([&](size_t, Index_ start, Index_ length) -> void {
            std::vector<Value_> vbuffer(ngenes);
            static_assert(!EigenMatrix_::IsRowMajor);
            auto optr = components.data() + static_cast<size_t>(start) * rank; // cast to avoid overflow.

            if (mat->is_sparse()) {
                std::vector<Index_> ibuffer(ngenes);
                auto ext = tatami::consecutive_extractor<true>(mat, false, start, length);

                for (Index_ c = start, end = start + length; c < end; ++c, optr += rank) {
                    auto range = ext->fetch(vbuffer.data(), ibuffer.data());
                    static_assert(!EigenMatrix_::IsRowMajor);
                    auto vptr = scaled_rotation.data();
                    for (size_t v = 0; v < rank; ++v, vptr += ngenes) {
                        auto& output = optr[v];
                        for (Index_ i = 0; i < range.number; ++i) {
                            output += vptr[range.index[i]] * range.value[i];
                        }
                    }
                }

            } else {
                auto ext = tatami::consecutive_extractor<false>(mat, false, start, length);
                for (Index_ c = start, end = start + length; c < end; ++c, optr += rank) {
                    auto ptr = ext->fetch(vbuffer.data()); 
                    static_assert(!EigenMatrix_::IsRowMajor);
                    auto vptr = scaled_rotation.data();
                    for (size_t v = 0; v < rank; ++v, vptr += ngenes) {
                        optr[v] = std::inner_product(vptr, vptr + ngenes, ptr, static_cast<Scalar>(0));
                    }
                }
            }
        }, ncells, nthreads);
    }
}

template<class EigenMatrix_, class EigenVector_>
void clean_up_projected(bool rows_are_dims, EigenMatrix_& projected, EigenVector_& D) {
    // Empirically centering to give nice centered PCs, because we can't
    // guarantee that the projection is centered in this manner.
    if (rows_are_dims) {
        for (size_t i = 0, iend = projected.rows(); i < iend; ++i) {
            projected.row(i).array() -= projected.row(i).sum() / projected.cols();
        }
    } else {
        for (size_t i = 0, iend = projected.cols(); i < iend; ++i) {
            projected.col(i).array() -= projected.col(i).sum() / projected.rows();
        }
    }

    // Just dividing by the number of observations - 1 regardless of weighting.
    typename EigenMatrix_::Scalar denom = (rows_are_dims ? projected.cols() : projected.rows()) - 1;
    for (auto& d : D) {
        d = d * d / denom;
    }
}

}

}

#endif
