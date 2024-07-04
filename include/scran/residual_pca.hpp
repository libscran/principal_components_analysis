#ifndef SCRAN_RESIDUAL_PCA_HPP
#define SCRAN_RESIDUAL_PCA_HPP

#include "tatami/tatami.hpp"
#include "irlba/irlba.hpp"
#include "Eigen/Dense"

#include <vector>
#include <cmath>

#include "pca_utils/general.hpp"
#include "pca_utils/blocking.hpp"

/**
 * @file residual_pca.hpp
 *
 * @brief Perform PCA on residuals after regressing out an uninteresting factor.
 */

namespace scran {

/**
 * @namespace scran::residual_pca
 * @brief Perform PCA on residuals after regressing out an uninteresting factor.
 *
 * A simple batch correction method involves centering the expression of each gene in each batch to remove systematic differences between batches.
 * The residuals are then used in PCA to obtain a batch-corrected low-dimensional representation of the dataset.
 * Unfortunately, naively centering the expression values will discard sparsity and reduce the computational efficiency of the PCA.
 * To avoid these drawbacks, `ResidualPca` defers the residual calculation until the matrix multiplication of the IRLBA step.
 * This yields the same results as the naive approach but is much faster as it can take advantage of efficient sparse operations.
 *
 * We can optionally scale each batch so that they contribute equally to the rotation vectors, regardless of their size.
 * This is achieved using the same approach described in `MultiBatchPca`,
 * whereby batches with more cells are downscaled during calculation of the rotation vectors.
 * The final PCs are then obtained by projecting the residuals onto the space defined by the rotation vectors.
 * This ensures that larger batches do not mask interesting variation in other batches.
 *
 * Note that the use of residuals for batch correction makes some strong assumptions about the batches,
 * e.g., they have the same composition and the batch effect is a consistent shift for all populations.
 * Non-linear correction algorithms are usually more effective, e.g., [MNN correction](https://github.com/LTLA/CppMnnCorrect).
 */
namespace residual_pca {

/**
 * @brief Options for `compute()`.
 */
struct Defaults {
    /**
     * Should genes be scaled to unit variance?
     * Genes with zero variance are ignored.
     */
    bool scale = false;

    /**
     * Should the PC matrix be transposed on output?
     * If `true`, the output matrix is column-major with cells in the columns, which is compatible with downstream **libscran** steps.
     */
    bool transpose = true;

    /**
     * Number of threads to use.
     */
    int num_threads = 1;

    /**
     * Policy to use for weighting batches of different size.
     */
    block_weights::Policy block_weight_policy = block_weights::Policy::VARIABLE;

    /**
     * Parameters for the variable block weights.
     * Only used when `Options::block_weight_policy = block_weights::Policy::VARIABLE`.
     */
    block_weights::VariableParameters variable_block_weight_parameters;

    /**
     * Whether to realize `tatami::Matrix` objects into an appropriate in-memory format before PCA.
     * This is typically faster but increases memory usage.
     */
    bool realize_matrix = true;

    /**
     * Further options to pass to `irlba::compute()`.
     */
    irlba::Options irlba_options;
};

/**
 * @cond
 */
namespace internal {

// This wrapper class mimics multiplication with the residuals,
// i.e., after subtracting the per-block mean from each cell.
template<class Matrix_, typename Block_, class EigenMatrix_, class EigenVector_>
struct RegressWrapper {
    RegressWrapper(const Matrix_& mat, const Block_* block, const EigenMatrix_& means) : my_mat(mat), my_block(block), my_means(means) {}

public:
    Eigen::Index rows() const { return my_mat.rows(); }
    Eigen::Index cols() const { return my_mat.cols(); }

public:
    struct Workspace {
        Workspace(size_t nblocks, irlba::WrappedWorkspace<Matrix_> c) : sub(nblocks), child(std::move(c)) {}
        EigenVector_ sub;
        EigenVector_ holding;
        irlba::WrappedWorkspace<Matrix_> child;
    };

    Workspace workspace() const {
        return Workspace(my_means.rows(), irlba::wrapped_workspace(my_mat));
    }

    template<class Right_>
    void multiply(const Right_& rhs, Workspace& work, EigenVector_& output) const {
        const auto& realized_rhs = [&]() {
            if constexpr(std::is_same<Right, EigenVector_>::value) {
                return rhs;
            } else {
                work.holding.noalias() = rhs;
                return work.holding;
            }
        }();

        irlba::wrapped_multiply(my_mat, realized_rhs, work.child, output);

        work.sub.noalias() = my_means * realized_rhs;
        for (Eigen::Index i = 0, end = output.size(); i < end; ++i) {
            auto& val = output.coeffRef(i);
            val -= work.sub.coeff(block[i]);
        }
    }

public:
    struct AdjointWorkspace {
        AdjointWorkspace(size_t nblocks, irlba::WrappedAdjointWorkspace<Matrix_> c) : aggr(nblocks), child(std::move(c)) {}
        EigenVector_ aggr;
        EigenVector_ holding;
        irlba::WrappedWorkspace<Matrix_> child;
    };

    AdjointWorkspace adjoint_workspace() const {
        return AdjointWorkspace(my_means.rows(), irlba::wrapped_adjoint_workspace(my_mat));
    }

    template<class Right_>
    void adjoint_multiply(const Right_& rhs, AdjointWorkspace& work, EigenVector_& output) const {
        const auto& realized_rhs = [&]() {
            if constexpr(std::is_same<Right, EigenVector_>::value) {
                return rhs;
            } else {
                work.holding.noalias() = rhs;
                return work.holding;
            }
        }();

        irlba::wrapped_adjoint_multiply(my_mat, realized_rhs, work.child, output);

        work.aggr.setZero();
        for (Eigen::Index i = 0, end = realized_rhs.size(); i < end; ++i) {
            work.aggr.coeffRef(my_block[i]) += realized_rhs.coeff(i); 
        }

        output.noalias() -= my_means.adjoint() * work.aggr;
    }

public:
    template<class EigenMatrix2_>
    EigenMatrix2_ realize() const {
        EigenMatrix2_ output = irlba::wrapped_realize<EigenMatrix2_>(mat);
        for (Eigen::Index r = 0, rend = output.rows(); r < rend; ++r) {
            output.row(r) -= my_means.row(block[r]);
        }
        return output;
    }

private:
    const Matrix_& mat;
    const Block_* block;
    const EigenMatrix_& means;
};

template<bool realize_matrix_, bool weight_, typename Value_, typename Index_, typename Block_, class EigenMatrix_, class EigenVector_>
void run_sparse(
    const tatami::Matrix<Value_, Index_>* mat, 
    const Block_* block, 
    const pca_utils::BlockingDetails& block_details, 
    int rank,
    const Options& options,
    EigenMatrix_& pcs, 
    EigenMatrix_& rotation, 
    EigenVector_& variance_explained, 
    EigenMatrix_& center_m,
    EigenVector_& scale_v,
    typename EigenVector_::Scalar& total_var) 
const {
    auto ngenes = mat->nrow(), ncells = mat->ncol(); 

    auto emat = [&]() {
        if constexpr(realize_matrix_) {
            // 'extracted' contains row-major contents...
            auto extracted = tatami::retrieve_compressed_sparse_contents<Value_, Index_>(
                mat, 
                /* row = */ true, 
                /* two_pass = */ false, 
                /* threads = */ options.num_threads
            );

            // But we effectively transpose it to CSC with genes in columns.
            Index_ ncells = mat->ncol();
            irlba::ParallelSparseMatrix emat(
                ncells,
                ngenes,
                std::move(extracted.value),
                std::move(extracted.index),
                std::move(extracted.pointers), 
                true,
                options.num_threads
            ); 
        } else {
            return pca_utils::TransposedTatamiWrapper<EigenVector_, Value_, Index_>(mat, options.num_threads), 
        }
    }();

    const auto& unwrapped = [&]() {
        if constexpr(realize_matrix_) {
            return emat;
        } else {
            return mat;
        }
    }();

    auto nblocks = block_details.num_blocks();
    center_m.resize(nblocks, ngenes);
    scale_v.resize(ngenes);
    pca_utils::compute_mean_and_variance_regress<weight_>(emat, block, block_details, center_m, scale_v, nthreads);
    total_var = pca_utils::process_scale_vector(scale, scale_v);

    RegressWrapper<decltype(emat), Block_, EigenMatrix_, EigenVector_> centered(emat, block, center_m);
    auto iopt = options.irlba_options;
    iopt.cap_number = true;

    if constexpr(weight) {
        if (scale) {
            irlba::Scaled<true, decltype(centered), EigenVector_> scaled(centered, scale_v, /* divide = */ true);
            irlba::Scaled<false, decltype(scaled), EigenVector_> weighted(scaled, block_details.expanded_weights, /* divide = */ false);
            irlba::compute(weighted, rank, components, rotation, variance_explained, iopt);
        } else {
            pca_utils::Scaled<false, decltype(centered), EigenVector_> weighted(&centered, block_details.expanded_weights, /* divide = */ false);
            irlba::compute(weighted, rank, components, rotation, variance_explained, iopt);
        }

        // This transposes 'pcs' to be a NDIM * NCELLS matrix.
        pca_utils::project_sparse_matrix(emat, pcs, rotation, scale, scale_v, nthreads);

        // Subtracting each block's mean from the PCs.
        EigenMatrix_ centering;
        if (scale) {
            centering = (center_m * (rotation.array().colwise() / scale_v.array()).matrix()).adjoint();
        } else {
            centering = (center_m * rotation).adjoint();
        }
        for (size_t i = 0, iend = pcs.cols(); i < iend; ++i) {
            pcs.col(i) -= centering.col(block[i]);
        }

        pca_utils::clean_up_projected<true>(pcs, variance_explained);
        if (!transpose) {
            pcs.adjointInPlace();
        }

    } else {
        if (scale) {
            irlba::Scaled<decltype(centered)> scaled(&centered, &scale_v);
            irb.run(scaled, pcs, rotation, variance_explained);
        } else {
            irb.run(centered, pcs, rotation, variance_explained);
        }

        pca_utils::clean_up(mat->ncol(), pcs, variance_explained);
        if (transpose) {
            pcs.adjointInPlace();
        }
    }
}

template<bool weight_, typename Data_, typename Index_, typename Block_>
void run_dense(
    const tatami::Matrix<Data_, Index_>* mat, 
    const Block_* block,
    const pca_utils::BlockingDetails<weight_>& block_details, 
    const irlba::Irlba& irb,
    EigenMatrix_& pcs, 
    EigenMatrix_& rotation, 
    EigenVector_& variance_explained, 
    EigenMatrix_& center_m,
    EigenVector_& scale_v,
    double& total_var) 
const {
    auto emat = pca_utils::extract_dense_for_pca(mat, nthreads); // get a column-major matrix with genes in columns.

    auto ngenes = emat.cols();
    auto nblocks = block_details.num_blocks();
    center_m.resize(nblocks, ngenes);
    scale_v.resize(ngenes);
    pca_utils::compute_mean_and_variance_regress<weight_>(emat, block, block_details, center_m, scale_v, nthreads);
    total_var = pca_utils::process_scale_vector(scale, scale_v);

    // Applying the centering and scaling directly so that we can run the PCA with no or fewer layers.
    tatami::parallelize([&](size_t, size_t start, size_t length) -> void {
        size_t ncells = emat.rows();
        double* ptr = emat.data() + static_cast<size_t>(start) * ncells;
        for (size_t g = start, end = start + length; g < end; ++g, ptr += ncells) {
            for (size_t c = 0; c < ncells; ++c) {
                ptr[c] -= center_m.coeff(block[c], g);
            }

            if (scale) {
                auto sd = scale_v[g];
                for (size_t c = 0; c < ncells; ++c) {
                    ptr[c] /= sd; // process_scale_vector should already protect against division by zero.
                }
            }
        }
    }, ngenes, nthreads);

    if constexpr(weight_) {
        pca_utils::SampleScaledWrapper<decltype(emat)> weighted(&emat, &(block_details.expanded_weights));
        irb.run(weighted, pcs, rotation, variance_explained);
        pcs.noalias() = emat * rotation;
        pca_utils::clean_up_projected<false>(pcs, variance_explained);
    } else {
        irb.run(emat, pcs, rotation, variance_explained);
        pca_utils::clean_up(pcs.rows(), pcs, variance_explained);
    }

    if (transpose) {
        pcs.adjointInPlace();
    }
}

template<typename Data_, typename Index_, typename Block_>
void dispatch(
    const tatami::Matrix<Data_, Index_>* mat, 
    const Block_* block,
    EigenMatrix_& pcs, 
    EigenMatrix_& rotation, 
    EigenVector_& variance_explained, 
    EigenMatrix_& center_m,
    EigenVector_& scale_v,
    double& total_var) 
const {
    irlba::EigenThreadScope t(nthreads);
    irlba::Irlba irb;
    irb.set_number(rank);
    irb.set_cap_number(true);

    if (block_weight_policy == WeightPolicy::NONE) {
        auto bdetails = pca_utils::compute_blocking_details(mat->ncol(), block);
        if (mat->sparse()) {
            run_sparse<false>(mat, block, bdetails, irb, pcs, rotation, variance_explained, center_m, scale_v, total_var);
        } else {
            run_dense<false>(mat, block, bdetails, irb, pcs, rotation, variance_explained, center_m, scale_v, total_var);
        }

    } else {
        auto bdetails = pca_utils::compute_blocking_details(mat->ncol(), block, block_weight_policy, variable_block_weight_parameters);
        if (mat->sparse()) {
            run_sparse<true>(mat, block, bdetails, irb, pcs, rotation, variance_explained, center_m, scale_v, total_var);
        } else {
            run_dense<true>(mat, block, bdetails, irb, pcs, rotation, variance_explained, center_m, scale_v, total_var);
        }
    }
}

}

/**
 * @brief Results of the PCA on the residuals.
 *
 * Instances should generally be constructed by the `compute()` function.
 */
template<typename EigenMatrix_, typename EigenVector_>
struct Results {
    /**
     * Matrix of principal components.
     * By default, each row corresponds to a PC while each column corresponds to a cell in the input matrix.
     * If `set_transpose()` is set to `false`, rows are cells instead.
     * The number of PCs is determined by `set_rank()`.
     */
    EigenMatrix_ components;

    /**
     * Variance explained by each PC.
     * Each entry corresponds to a column in `pcs` and is in decreasing order.
     */
    EigenVector_ variance_explained;

    /**
     * Total variance of the dataset (possibly after scaling, if `set_scale()` is set to `true`).
     * This can be used to divide `variance_explained` to obtain the percentage of variance explained.
     */
    typename EigenVector_::Scalar total_variance = 0;

    /**
     * Rotation matrix, only returned if `ResidualPca::set_return_rotation()` is `true`.
     * Each row corresponds to a feature while each column corresponds to a PC.
     * The number of PCs is determined by `set_rank()`.
     * If feature filtering was performed, the number of rows is equal to the number of features remaining after filtering.
     */
    EigenMatrix_ rotation;

    /**
     * Centering matrix, only returned if `ResidualPca::set_return_center()` is `true`.
     * Each row corresponds to a row in the input matrix and each column corresponds to a block, 
     * such that each entry contains the mean for a particular feature in the corresponding block.
     * If feature filtering was performed, the number of rows is equal to the number of features remaining after filtering.
     */
    EigenMatrix_ center;

    /**
     * Scaling vector, only returned if `ResidualPca::set_return_center()` is `true`.
     * Each entry corresponds to a row in the input matrix and contains the scaling factor used to divide the feature values if `ResidualPca::set_scale()` is `true`.
     * If feature filtering was performed, the length is equal to the number of features remaining after filtering.
     */
    EigenVector_ scale;
};

/**
 * Run the blocked PCA on an input gene-by-cell matrix.
 *
 * @tparam Data_ Floating point type for the data.
 * @tparam Index_ Integer type for the indices.
 * @tparam Block_ Integer type for the blocking factor.
 *
 * @param[in] mat Pointer to the input matrix.
 * Columns should contain cells while rows should contain genes.
 * @param[in] block Pointer to an array of length equal to the number of cells, 
 * containing the block assignment for each cell - see `count_blocks()` for details.
 *
 * @return A `Results` object containing the PCs and the variance explained.
 */
template<typename Data_, typename Index_, typename Block_>
Results compute(const tatami::Matrix<Data_, Index_>* mat, const Block_* block) const {
    Results output;
    internal::dispatch(mat, block, output.components, output.rotation, output.variance_explained, output.center, output.scale, output.total_variance);
    return output;
}

}

#endif
