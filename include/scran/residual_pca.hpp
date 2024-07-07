#ifndef SCRAN_RESIDUAL_PCA_HPP
#define SCRAN_RESIDUAL_PCA_HPP

#include "tatami/tatami.hpp"
#include "irlba/irlba.hpp"
#include "irlba/parallel.hpp"
#include "Eigen/Dense"

#include <vector>
#include <cmath>

#include "pca_utils/general.hpp"
#include "pca_utils/blocking.hpp"
#include "pca_utils/TransposedTatamiWrapper.hpp"

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
struct Options {
    /**
     * @cond
     */
    Options() {
        irlba_options.cap_number = true;
    }
    /**
     * @endcond
     */

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
class RegressWrapper {
public:
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
            if constexpr(std::is_same<Right_, EigenVector_>::value) {
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
            val -= work.sub.coeff(my_block[i]);
        }
    }

public:
    struct AdjointWorkspace {
        AdjointWorkspace(size_t nblocks, irlba::WrappedAdjointWorkspace<Matrix_> c) : aggr(nblocks), child(std::move(c)) {}
        EigenVector_ aggr;
        EigenVector_ holding;
        irlba::WrappedAdjointWorkspace<Matrix_> child;
    };

    AdjointWorkspace adjoint_workspace() const {
        return AdjointWorkspace(my_means.rows(), irlba::wrapped_adjoint_workspace(my_mat));
    }

    template<class Right_>
    void adjoint_multiply(const Right_& rhs, AdjointWorkspace& work, EigenVector_& output) const {
        const auto& realized_rhs = [&]() {
            if constexpr(std::is_same<Right_, EigenVector_>::value) {
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
        EigenMatrix2_ output = irlba::wrapped_realize<EigenMatrix2_>(my_mat);
        for (Eigen::Index r = 0, rend = output.rows(); r < rend; ++r) {
            output.row(r) -= my_means.row(my_block[r]);
        }
        return output;
    }

private:
    const Matrix_& my_mat;
    const Block_* my_block;
    const EigenMatrix_& my_means;
};

template<typename Block_, class EigenMatrix_>
void subtract_centers_from_components(const Block_* block, EigenMatrix_& components, const EigenMatrix_& center_m, const EigenMatrix_& rotation) {
    EigenMatrix_ centering = (center_m * rotation).adjoint();
    for (size_t i = 0, iend = components.cols(); i < iend; ++i) {
        components.col(i) -= centering.col(block[i]);
    }
}

template<bool realize_matrix_, typename Value_, typename Index_, typename Block_, class EigenMatrix_, class EigenVector_>
void run_sparse(
    const tatami::Matrix<Value_, Index_>* mat, 
    const Block_* block, 
    const pca_utils::BlockingDetails<Index_, EigenVector_>& block_details, 
    int rank,
    const Options& options,
    EigenMatrix_& components, 
    EigenMatrix_& rotation, 
    EigenVector_& variance_explained, 
    EigenMatrix_& center_m,
    EigenVector_& scale_v,
    typename EigenVector_::Scalar& total_var) 
{
    Index_ ngenes = mat->nrow(), ncells = mat->ncol(); 

    auto emat = [&]{
        if constexpr(realize_matrix_) {
            // 'extracted' contains row-major contents...
            auto extracted = tatami::retrieve_compressed_sparse_contents<Value_, Index_>(
                mat, 
                /* row = */ true, 
                /* two_pass = */ false, 
                /* threads = */ options.num_threads
            );

            // But we effectively transpose it to CSC with genes in columns.
            return irlba::ParallelSparseMatrix(
                ncells,
                ngenes,
                std::move(extracted.value),
                std::move(extracted.index),
                std::move(extracted.pointers), 
                true,
                options.num_threads
            ); 
        } else {
            return pca_utils::TransposedTatamiWrapper<EigenVector_, Value_, Index_>(mat, options.num_threads);
        }
    }();

    auto nblocks = block_details.block_size.size();
    center_m.resize(nblocks, ngenes);
    scale_v.resize(ngenes);
    if constexpr(realize_matrix_) {
        pca_utils::compute_blockwise_mean_and_variance_realized_sparse(emat, block, block_details, center_m, scale_v, options.num_threads);
    } else {
        pca_utils::compute_blockwise_mean_and_variance_tatami(mat, block, block_details, center_m, scale_v, options.num_threads);
    }
    total_var = pca_utils::process_scale_vector(options.scale, scale_v);

    RegressWrapper<decltype(emat), Block_, EigenMatrix_, EigenVector_> centered(emat, block, center_m);

    if (block_details.weighted) {
        if (options.scale) {
            irlba::Scaled<true, decltype(centered), EigenVector_> scaled(centered, scale_v, /* divide = */ true);
            irlba::Scaled<false, decltype(scaled), EigenVector_> weighted(scaled, block_details.expanded_weights, /* divide = */ false);
            irlba::compute(weighted, rank, components, rotation, variance_explained, options.irlba_options);
        } else {
            irlba::Scaled<false, decltype(centered), EigenVector_> weighted(centered, block_details.expanded_weights, /* divide = */ false);
            irlba::compute(weighted, rank, components, rotation, variance_explained, options.irlba_options);
        }

        EigenMatrix_ tmp;
        const auto& scaled_rotation = pca_utils::scale_rotation_matrix(rotation, options.scale, scale_v, tmp);

        // This transposes 'components' to be a NDIM * NCELLS matrix.
        if constexpr(realize_matrix_) {
            pca_utils::project_matrix_realized_sparse(emat, components, scaled_rotation, options.num_threads);
        } else {
            pca_utils::project_matrix_transposed_tatami(mat, components, scaled_rotation, options.num_threads);
        }

        // Subtracting each block's mean from the PCs.
        subtract_centers_from_components(block, components, center_m, scaled_rotation);

        pca_utils::clean_up_projected(/* rows_are_dims = */ true, components, variance_explained);
        if (!options.transpose) {
            components.adjointInPlace();
        }

    } else {
        if (options.scale) {
            irlba::Scaled<true, decltype(centered), EigenVector_> scaled(centered, scale_v, /* divide = */ true);
            irlba::compute(scaled, rank, components, rotation, variance_explained, options.irlba_options);
        } else {
            irlba::compute(centered, rank, components, rotation, variance_explained, options.irlba_options);
        }

        pca_utils::clean_up(mat->ncol(), components, variance_explained);
        if (options.transpose) {
            components.adjointInPlace();
        }
    }
}

template<bool realize_matrix_, typename Value_, typename Index_, typename Block_, class EigenMatrix_, class EigenVector_>
void run_dense(
    const tatami::Matrix<Value_, Index_>* mat, 
    const Block_* block,
    const pca_utils::BlockingDetails<Index_, EigenVector_>& block_details, 
    int rank,
    const Options& options,
    EigenMatrix_& components, 
    EigenMatrix_& rotation, 
    EigenVector_& variance_explained, 
    EigenMatrix_& center_m,
    EigenVector_& scale_v,
    double& total_var) 
{
    Index_ ngenes = mat->nrow(), ncells = mat->ncol(); 

    auto emat = [&]() {
        if constexpr(realize_matrix_) {
            // Perform an implicit transposition by performing a row-major extraction
            // into a column-major transposed matrix.
            EigenMatrix_ emat(ncells, ngenes); 
            static_assert(!EigenMatrix_::IsRowMajor);
            tatami::convert_to_dense(mat, /* row_major = */ true, emat.data(), options.num_threads);
            return emat;
        } else {
            return pca_utils::TransposedTatamiWrapper<EigenVector_, Value_, Index_>(mat, options.num_threads);
        }
    }();

    auto nblocks = block_details.block_size.size();
    center_m.resize(nblocks, ngenes);
    scale_v.resize(ngenes);
    if constexpr(realize_matrix_) {
        pca_utils::compute_blockwise_mean_and_variance_realized_dense(emat, block, block_details, center_m, scale_v, options.num_threads);
    } else {
        pca_utils::compute_blockwise_mean_and_variance_tatami(mat, block, block_details, center_m, scale_v, options.num_threads);
    }
    total_var = pca_utils::process_scale_vector(options.scale, scale_v);

    if constexpr(realize_matrix_) {
        // Applying the centering and scaling directly so that we can run the PCA with no or fewer layers.
        for (Index_ c = 0; c < ncells; ++c) {
            emat.row(c) -= center_m.row(block[c]);
        }
        if (options.scale) {
            emat.array().rowwise() /= scale_v.adjoint().array(); // process_scale_vector should already protect against division by zero.
        }

        if (block_details.weighted) {
            irlba::Scaled<false, decltype(emat), EigenVector_> weighted(emat, block_details.expanded_weights, /* divide = */ false);
            irlba::compute(weighted, rank, components, rotation, variance_explained, options.irlba_options);
            components.noalias() = emat * rotation;
            pca_utils::clean_up_projected(/* rows_are_dims = */ false, components, variance_explained);
        } else {
            irlba::compute(emat, rank, components, rotation, variance_explained, options.irlba_options);
            pca_utils::clean_up(components.rows(), components, variance_explained);
        }

        if (options.transpose) {
            components.adjointInPlace();
        }

    } else {
        RegressWrapper<decltype(emat), Block_, EigenMatrix_, EigenVector_> centered(emat, block, center_m);

        if (block_details.weighted) {
            if (options.scale) {
                irlba::Scaled<true, decltype(centered), EigenVector_> scaled(centered, scale_v, /* divide = */ true);
                irlba::Scaled<false, decltype(scaled), EigenVector_> weighted(scaled, block_details.expanded_weights, /* divide = */ false);
                irlba::compute(weighted, rank, components, rotation, variance_explained, options.irlba_options);
            } else {
                irlba::Scaled<false, decltype(centered), EigenVector_> weighted(centered, block_details.expanded_weights, /* divide = */ false);
                irlba::compute(weighted, rank, components, rotation, variance_explained, options.irlba_options);
            }

            EigenMatrix_ tmp;
            const auto& scaled_rotation = pca_utils::scale_rotation_matrix(rotation, options.scale, scale_v, tmp);

            // This transposes 'components' to be a NDIM * NCELLS matrix.
            if constexpr(realize_matrix_) {
                components.noalias() = (emat * scaled_rotation).adjoint();
            } else {
                pca_utils::project_matrix_transposed_tatami(mat, components, scaled_rotation, options.num_threads);
            }

            // Subtracting each block's mean from the PCs.
            subtract_centers_from_components(block, components, center_m, scaled_rotation);

            pca_utils::clean_up_projected(/* rows_are_dims = */ true, components, variance_explained);
            if (!options.transpose) {
                components.adjointInPlace();
            }

        } else {
            if (options.scale) {
                irlba::Scaled<true, decltype(centered), EigenVector_> scaled(centered, scale_v, /* divide = */ true);
                irlba::compute(scaled, rank, components, rotation, variance_explained, options.irlba_options);
            } else {
                irlba::compute(centered, rank, components, rotation, variance_explained, options.irlba_options);
            }

            pca_utils::clean_up(components.rows(), components, variance_explained);
            if (options.transpose) {
                components.adjointInPlace();
            }
        }
    }
}

template<typename Value_, typename Index_, typename Block_, class EigenMatrix_, class EigenVector_>
void dispatch(
    const tatami::Matrix<Value_, Index_>* mat, 
    const Block_* block,
    int rank,
    const Options& options,
    EigenMatrix_& components, 
    EigenMatrix_& rotation, 
    EigenVector_& variance_explained, 
    EigenMatrix_& center_m,
    EigenVector_& scale_v,
    double& total_var) 
{
    irlba::EigenThreadScope t(options.num_threads);
    auto bdetails = pca_utils::compute_blocking_details<EigenVector_>(mat->ncol(), block, options.block_weight_policy, options.variable_block_weight_parameters);

    if (mat->sparse()) {
        if (options.realize_matrix) {
            run_sparse<true>(mat, block, bdetails, rank, options, components, rotation, variance_explained, center_m, scale_v, total_var);
        } else {
            run_sparse<false>(mat, block, bdetails, rank, options, components, rotation, variance_explained, center_m, scale_v, total_var);
        }
    } else {
        if (options.realize_matrix) {
            run_dense<true>(mat, block, bdetails, rank, options, components, rotation, variance_explained, center_m, scale_v, total_var);
        } else {
            run_dense<false>(mat, block, bdetails, rank, options, components, rotation, variance_explained, center_m, scale_v, total_var);
        }
    }
}

}
/**
 * @endcond
 */

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
     * If `Options::transpose = false`, rows are cells instead.
     * The number of PCs is determined by the `rank` used in `compute()`.
     */
    EigenMatrix_ components;

    /**
     * Variance explained by each PC.
     * Each entry corresponds to a column in `components` and is in decreasing order.
     */
    EigenVector_ variance_explained;

    /**
     * Total variance of the dataset (possibly after scaling, if `Options::scale = true`).
     * This can be used to divide `variance_explained` to obtain the percentage of variance explained.
     */
    typename EigenVector_::Scalar total_variance = 0;

    /**
     * Rotation matrix.
     * Each row corresponds to a gene while each column corresponds to a PC.
     * The number of PCs is determined by the `rank` used in `compute()`.
     */
    EigenMatrix_ rotation;

    /**
     * Centering matrix.
     * Each row corresponds to a block and each column corresponds to a gene.
     * Each entry contains the mean for a particular gene in the corresponding block.
     */
    EigenMatrix_ center;

    /**
     * Scaling vector, only returned if `Options::scale = true`.
     * Each entry corresponds to a row in the input matrix and contains the scaling factor used to divide that gene's values if `Options::scale = true`.
     */
    EigenVector_ scale;
};

/**
 * Run PCA on the residuals after regressing out a blocking factor from the rows of a gene-by-cell matrix.
 *
 * @tparam EigenMatrix_ A floating-point `Eigen::Matrix` class.
 * @tparam EigenVector_ A floating-point `Eigen::Vector` class.
 * @tparam Value_ Type of the matrix data.
 * @tparam Index_ Integer type for the indices.
 * @tparam Block_ Integer type for the blocking factor.
 *
 * @param[in] mat Pointer to the input matrix.
 * Columns should contain cells while rows should contain genes.
 * @param[in] block Pointer to an array of length equal to the number of cells, 
 * containing the block assignment for each cell. 
 * Each assignment should be an integer in \f$[0, N)\f$ where \f$N\f$ is the number of blocks.
 * @param rank Number of PCs to compute.
 * This should be no greater than the maximum number of PCs, i.e., the smaller dimension of the input matrix;
 * otherwise, only the maximum number of PCs will be reported in the results.
 * @param options Further options.
 *
 * @return The results of the PCA on the residuals. 
 */
template<typename EigenMatrix_ = Eigen::MatrixXd, class EigenVector_ = Eigen::VectorXd, typename Value_ = double, typename Index_ = int, typename Block_ = int>
Results<EigenMatrix_, EigenVector_> compute(const tatami::Matrix<Value_, Index_>* mat, const Block_* block, int rank, const Options& options) {
    Results<EigenMatrix_, EigenVector_> output;
    internal::dispatch(mat, block, rank, options, output.components, output.rotation, output.variance_explained, output.center, output.scale, output.total_variance);
    return output;
}

}

}

#endif
