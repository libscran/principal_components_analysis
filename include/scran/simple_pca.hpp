#ifndef SCRAN_SIMPLE_PCA_HPP
#define SCRAN_SIMPLE_PCA_HPP

#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"
#include "irlba/irlba.hpp"
#include "irlba/parallel.hpp"
#include "Eigen/Dense"

#include <vector>
#include <cmath>

#include "pca_utils.hpp"

/**
 * @file simple_pca.hpp
 * @brief Perform a simple PCA on a gene-by-cell matrix.
 */

namespace scran {

/**
 * @namespace scran::simple_pca
 * @brief Perform a simple PCA on a gene-cell matrix.
 *
 * Principal components analysis (PCA) is a helpful technique for data compression and denoising.
 * The idea is that the earlier PCs capture most of the systematic biological variation while the later PCs capture random technical noise.
 * Thus, we can reduce the size of the data and eliminate noise by only using the earlier PCs for further analyses.
 * Most practitioners will keep the first 10-50 PCs, though the exact choice is fairly arbitrary.
 */
namespace simple_pca {

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

template<bool sparse_, typename Value_, typename Index_, class EigenVector_>
void compute_row_means_and_variances(const tatami::Matrix<Value_, Index_>* mat, int num_threads, EigenVector_& center_v, EigenVector_& scale_v) {
    if (mat->prefer_rows()) {
        tatami::parallelize([&](size_t, Index_ start, Index_ length) -> void {
            tatami::Options opt;
            opt.sparse_extract_index = false;
            auto ext = tatami::consecutive_extractor<sparse_>(mat, true, start, length, opt);
            auto ncells = mat->ncol();
            std::vector<Value_> vbuffer(ncells);

            for (Index_ r = start, end = start + length; r < end; ++r) {
                auto results = [&]() {
                    if constexpr(sparse_) {
                        auto range = ext->fetch(vbuffer.data(), NULL);
                        return tatami_stats::variances::direct(range.value, range.number, ncells, /* skip_nan = */ false);
                    } else {
                        auto ptr = ext->fetch(vbuffer.data());
                        return tatami_stats::variances::direct(ptr, ncells, /* skip_nan = */ false);
                    }
                }();
                center_v.coeffRef(r) = results.first;
                scale_v.coeffRef(r) = results.second;
            }
        }, mat->nrow(), num_threads);

    } else {
        tatami::parallelize([&](size_t t, Index_ start, Index_ length) -> void {
            tatami::Options opt;
            auto ncells = mat->ncol();
            auto ext = tatami::consecutive_extractor<sparse_>(mat, false, 0, ncells, start, length, opt);

            typedef typename EigenVector_::Scalar Scalar;
            tatami_stats::LocalOutputBuffer<Scalar> cbuffer(t, start, length, center_v.data());
            tatami_stats::LocalOutputBuffer<Scalar> sbuffer(t, start, length, scale_v.data());

            auto running = [&]() {
                if constexpr(sparse_) {
                    return tatami_stats::variances::RunningSparse<Scalar, Value_, Index_>(length, cbuffer.data(), sbuffer.data(), /* skip_nan = */ false, /* subtract = */ start);
                } else {
                    return tatami_stats::variances::RunningDense<Scalar, Value_, Index_>(length, cbuffer.data(), sbuffer.data(), /* skip_nan = */ false);
                }
            }();

            std::vector<Value_> vbuffer(length);
            typename std::conditional<sparse_, std::vector<Index_>, Index_>::type ibuffer(length);
            for (Index_ r = 0; r < ncells; ++r) {
                if constexpr(sparse_) {
                    auto range = ext->fetch(vbuffer.data(), ibuffer.data());
                    running.add(range.value, range.index, range.number);
                } else {
                    auto ptr = ext->fetch(vbuffer.data());
                    running.add(ptr);
                }
            }

            running.finish();
            cbuffer.transfer();
            sbuffer.transfer();
        }, mat->nrow(), num_threads);
    }
}

template<class IrlbaMatrix_, class EigenMatrix_, class EigenVector_>
void run_irlba_deferred(
    const IrlbaMatrix_& mat,
    int rank,
    const Options& options,
    EigenMatrix_& components, 
    EigenMatrix_& rotation, 
    EigenVector_& variance_explained,
    EigenVector_& center_v,
    EigenVector_& scale_v)
{
    irlba::Centered<IrlbaMatrix_, EigenVector_> centered(mat, center_v);
    if (options.scale) {
        irlba::Scaled<true, decltype(centered), EigenVector_> scaled(centered, scale_v, true);
        irlba::compute(scaled, rank, components, rotation, variance_explained, options.irlba_options);
    } else {
        irlba::compute(centered, rank, components, rotation, variance_explained, options.irlba_options);
    }
}

template<typename Value_, typename Index_, class EigenMatrix_, class EigenVector_>
void run_sparse(
    const tatami::Matrix<Value_, Index_>* mat, 
    int rank,
    const Options& options,
    EigenMatrix_& components, 
    EigenMatrix_& rotation, 
    EigenVector_& variance_explained,
    EigenVector_& center_v,
    EigenVector_& scale_v,
    typename EigenVector_::Scalar& total_var)
{
    Index_ ngenes = mat->nrow();
    center_v.resize(ngenes);
    scale_v.resize(ngenes);

    if (options.realize_matrix) {
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

        tatami::parallelize([&](size_t, size_t start, size_t length) -> void {
            const auto& ptrs = emat.get_pointers();
            const auto& values = emat.get_values();
            for (size_t r = start, end = start + length; r < end; ++r) {
                auto offset = ptrs[r];
                Index_ num_nonzero = ptrs[r + 1] - offset;
                auto results = tatami_stats::variances::direct(values.data() + offset, num_nonzero, ncells, /* skip_nan = */ false);
                center_v.coeffRef(r) = results.first;
                scale_v.coeffRef(r) = results.second;
            }
        }, ngenes, options.num_threads);

        total_var = pca_utils::process_scale_vector(options.scale, scale_v);
        run_irlba_deferred(emat, rank, options, components, rotation, variance_explained, center_v, scale_v);

    } else {
        compute_row_means_and_variances<true>(mat, options.num_threads, center_v, scale_v);
        total_var = pca_utils::process_scale_vector(options.scale, scale_v);
        run_irlba_deferred(
            pca_utils::TransposedTatamiWrapper<EigenVector_, Value_, Index_>(mat, options.num_threads), 
            rank, 
            options, 
            components, 
            rotation, 
            variance_explained, 
            center_v, 
            scale_v
        );
    }
}

template<typename Value_, typename Index_, class EigenMatrix_, class EigenVector_>
void run_dense(
    const tatami::Matrix<Value_, Index_>* mat, 
    int rank,
    const Options& options,
    EigenMatrix_& components, 
    EigenMatrix_& rotation, 
    EigenVector_& variance_explained, 
    EigenVector_& center_v,
    EigenVector_& scale_v,
    typename EigenVector_::Scalar& total_var)
{
    Index_ ngenes = mat->nrow();
    center_v.resize(ngenes);
    scale_v.resize(ngenes);

    if (options.realize_matrix) {
        // Create a matrix with genes in columns.
        Index_ ncells = mat->ncol();
        EigenMatrix_ emat(ncells, ngenes);

        // If emat is row-major, we want to fill it with columns of 'mat', so row_major = false.
        // If emat is column-major, we want to fill it with rows of 'mat', so row_major = true.
        tatami::convert_to_dense(mat, /* row_major = */ !emat.IsRowMajor, emat.data(), options.num_threads);

        center_v.array() = emat.array().colwise().sum();
        if (ncells) {
            center_v /= ncells;
        } else {
            std::fill(center_v.begin(), center_v.end(), std::numeric_limits<typename EigenVector_::Scalar>::quiet_NaN());
        }
        emat.array().rowwise() -= center_v.adjoint().array(); // applying it to avoid wasting time with deferred operations inside IRLBA.

        scale_v.array() = emat.array().colwise().squaredNorm();
        if (ncells > 1) {
            scale_v /= ncells - 1;
        } else {
            std::fill(scale_v.begin(), scale_v.end(), std::numeric_limits<typename EigenVector_::Scalar>::quiet_NaN());
        }

        total_var = pca_utils::process_scale_vector(options.scale, scale_v);
        if (options.scale) {
            emat.array().rowwise() /= scale_v.adjoint().array();
        }

        irlba::compute(emat, rank, components, rotation, variance_explained, options.irlba_options);

    } else {
        compute_row_means_and_variances<false>(mat, options.num_threads, center_v, scale_v);
        total_var = pca_utils::process_scale_vector(options.scale, scale_v);
        run_irlba_deferred(
            pca_utils::TransposedTatamiWrapper<EigenVector_, Value_, Index_>(mat, options.num_threads), 
            rank, 
            options, 
            components, 
            rotation, 
            variance_explained, 
            center_v, 
            scale_v
        );
    }
}

template<typename Value_, typename Index_, class EigenMatrix_, class EigenVector_>
void dispatch(
    const tatami::Matrix<Value_, Index_>* mat, 
    int rank, 
    const Options& options,
    EigenMatrix_& components, 
    EigenMatrix_& rotation, 
    EigenVector_& variance_explained, 
    EigenVector_& center_v,
    EigenVector_& scale_v,
    typename EigenVector_::Scalar& total_var) 
{
    irlba::EigenThreadScope t(options.num_threads);

    if (mat->sparse()) {
        run_sparse(mat, rank, options, components, rotation, variance_explained, center_v, scale_v, total_var);
    } else {
        run_dense(mat, rank, options, components, rotation, variance_explained, center_v, scale_v, total_var);
    }

    pca_utils::clean_up(mat->ncol(), components, variance_explained);
    if (options.transpose) {
        components.adjointInPlace();
    }
}

}
/**
 * @cond
 */

/**
 * @brief Container for the PCA results.
 * @tparam EigenMatrix_ A floating-point `Eigen::Matrix` class.
 * @tparam EigenVector_ A floating-point `Eigen::Vector` class.
 *
 * Instances should be constructed by the `compute()` function.
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
     * Each row corresponds to a feature while each column corresponds to a PC.
     * The number of PCs is determined by the `rank` used in `compute()`.
     */
    EigenMatrix_ rotation;

    /**
     * Centering vector.
     * Each entry corresponds to a row in the matrix and contains the mean value for that feature.
     */
    EigenVector_ center;

    /**
     * Scaling vector, only returned if `Options::scale = true`.
     * Each entry corresponds to a row in the matrix and contains the scaling factor used to divide the feature values if `Options::scale = true`.
     */
    EigenVector_ scale;
};

/**
 * Run PCA on an input gene-by-cell matrix.
 *
 * @tparam EigenMatrix_ A floating-point `Eigen::Matrix` class.
 * @tparam EigenVector_ A floating-point `Eigen::Vector` class.
 * @tparam Value_ Type of the matrix data.
 * @tparam Index_ Integer type for the indices.
 *
 * @param[in] mat Pointer to the input matrix.
 * Columns should contain cells while rows should contain genes.
 * @param rank Number of PCs to compute.
 * This should be no greater than the maximum number of PCs, i.e., the smaller dimension of the input matrix;
 * otherwise, only the maximum number of PCs will be reported in the results.
 * @param options Further options.
 *
 * @return The results of the PCA on `mat`.
 */
template<typename EigenMatrix_ = Eigen::MatrixXd, class EigenVector_ = Eigen::VectorXd, typename Value_ = double, typename Index_ = int>
Results<EigenMatrix_, EigenVector_> compute(const tatami::Matrix<Value_, Index_>* mat, int rank, const Options& options) {
    Results<EigenMatrix_, EigenVector_> output;
    internal::dispatch(mat, rank, options, output.components, output.rotation, output.variance_explained, output.center, output.scale, output.total_variance);

    if (!options.scale) {
        output.scale = EigenVector_();
    }

    return output;
}

}

}

#endif
