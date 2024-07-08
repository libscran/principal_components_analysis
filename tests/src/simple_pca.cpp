#include <gtest/gtest.h>

#include "simulate_vector.h"
#include "compare_pcs.h"

#include "tatami/tatami.hpp"

#include "simple_pca.hpp"

class SimplePcaTestCore {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_row, dense_column, sparse_row, sparse_column;

    static void assemble() {
        if (dense_row) {
            return;
        }

        size_t nr = 199, nc = 165;
        auto vec = simulate_sparse_vector(nr * nc, 0.1, /* lower = */ -10, /* upper = */ 10, /* seed = */ 69);
        dense_row.reset(new tatami::DenseRowMatrix<double, int>(nr, nc, std::move(vec)));
        dense_column = tatami::convert_to_dense(dense_row.get(), false);
        sparse_row = tatami::convert_to_compressed_sparse(dense_row.get(), true);
        sparse_column = tatami::convert_to_compressed_sparse(dense_row.get(), false);
    }
};

/******************************************/

class SimplePcaBasicTest : public ::testing::TestWithParam<std::tuple<bool, int, int> >, public SimplePcaTestCore {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(SimplePcaBasicTest, Test) {
    auto param = GetParam();
    bool scale = std::get<0>(param);
    int rank = std::get<1>(param);
    int threads = std::get<2>(param);

    scran::simple_pca::Options opt;
    opt.scale = scale;
    auto ref = scran::simple_pca::compute(dense_row.get(), rank, opt);

    if (threads == 1) {
        EXPECT_EQ(ref.variance_explained.size(), rank);
        EXPECT_EQ(ref.components.rows(), rank);
        EXPECT_EQ(ref.components.cols(), dense_row->ncol());

        // Checking that we scaled the PCs correctly.
        are_pcs_centered(ref.components);
        size_t NC = dense_row->ncol();
        for (int r = 0; r < rank; ++r) {
            double var = 0;
            auto ptr = ref.components.data() + r;
            for (size_t c = 0; c < NC; ++c, ptr += rank) {
                var += (*ptr) * (*ptr);
            }
            var /= NC - 1;
            EXPECT_FLOAT_EQ(var, ref.variance_explained[r]);
        }

        if (scale) {
            EXPECT_FLOAT_EQ(dense_row->nrow(), ref.total_variance);
        } else {
            auto vars = tatami_stats::variances::by_row(dense_row.get());
            auto total_var = std::accumulate(vars.begin(), vars.end(), 0.0);
            EXPECT_FLOAT_EQ(total_var, ref.total_variance);
        }

        EXPECT_TRUE(ref.total_variance >= std::accumulate(ref.variance_explained.begin(), ref.variance_explained.end(), 0.0));

        if (scale) {
            EXPECT_EQ(ref.scale.size(), dense_row->nrow());
        } else {
            EXPECT_EQ(ref.scale.size(), 0);
        }

    } else {
        // Results should be EXACTLY the same with parallelization.
        opt.num_threads = threads;
        auto res1 = scran::simple_pca::compute(dense_row.get(), rank, opt);
        EXPECT_EQ(ref.components, res1.components);
        EXPECT_EQ(ref.variance_explained, res1.variance_explained);
        EXPECT_EQ(ref.total_variance, res1.total_variance);
    }

    // Checking that we get more-or-less the same results. 
    auto res2 = scran::simple_pca::compute(dense_column.get(), rank, opt);
    expect_equal_pcs(ref.components, res2.components);
    expect_equal_vectors(ref.variance_explained, res2.variance_explained);
    EXPECT_FLOAT_EQ(ref.total_variance, res2.total_variance);

    auto res3 = scran::simple_pca::compute(sparse_row.get(), rank, opt);
    expect_equal_pcs(ref.components, res3.components);
    expect_equal_vectors(ref.variance_explained, res3.variance_explained);
    EXPECT_FLOAT_EQ(ref.total_variance, res3.total_variance);

    auto res4 = scran::simple_pca::compute(sparse_column.get(), rank, opt);
    expect_equal_pcs(ref.components, res4.components);
    expect_equal_vectors(ref.variance_explained, res4.variance_explained);
    EXPECT_FLOAT_EQ(ref.total_variance, res4.total_variance);

    // Checking that we get more-or-less the same results. 
    opt.realize_matrix = false;
    auto tres1 = scran::simple_pca::compute(dense_row.get(), rank, opt);
    expect_equal_pcs(ref.components, tres1.components);
    expect_equal_vectors(ref.variance_explained, tres1.variance_explained);
    EXPECT_FLOAT_EQ(ref.total_variance, tres1.total_variance);

    auto tres2 = scran::simple_pca::compute(dense_column.get(), rank, opt);
    expect_equal_pcs(ref.components, tres2.components);
    expect_equal_vectors(ref.variance_explained, tres2.variance_explained);
    EXPECT_FLOAT_EQ(ref.total_variance, tres2.total_variance);

    auto tres3 = scran::simple_pca::compute(sparse_row.get(), rank, opt);
    expect_equal_pcs(ref.components, tres3.components);
    expect_equal_vectors(ref.variance_explained, tres3.variance_explained);
    EXPECT_FLOAT_EQ(ref.total_variance, tres3.total_variance);

    auto tres4 = scran::simple_pca::compute(sparse_column.get(), rank, opt);
    expect_equal_pcs(ref.components, tres4.components);
    expect_equal_vectors(ref.variance_explained, tres4.variance_explained);
    EXPECT_FLOAT_EQ(ref.total_variance, tres4.total_variance);
}

INSTANTIATE_TEST_SUITE_P(
    SimplePca,
    SimplePcaBasicTest,
    ::testing::Combine(
        ::testing::Values(false, true), // to scale or not to scale?
        ::testing::Values(2, 5, 10), // number of PCs to obtain
        ::testing::Values(1, 3) // number of threads
    )
);

/******************************************/

class SimplePcaMoreTest : public ::testing::TestWithParam<std::tuple<bool, int> >, public SimplePcaTestCore {};

TEST_P(SimplePcaMoreTest, ZeroVariance) {
    auto param = GetParam();
    bool scale = std::get<0>(param);
    int rank = std::get<1>(param);

    size_t nr = 109, nc = 153;
    auto vec = simulate_sparse_vector(nr * nc, 0.1, /* lower = */ -10, /* upper = */ 10, /* seed = */ scale + rank + 100);

    auto copy = vec;
    size_t last_row = (nr - 1) * nc;
    std::fill(copy.begin() + last_row, copy.begin() + last_row + nc, 0);
    tatami::DenseRowMatrix<double, int> has_zero(nr, nc, std::move(copy));

    std::vector<double> removed(vec.begin(), vec.begin() + last_row);
    tatami::DenseRowMatrix<double, int> leftovers(nr - 1, nc, std::move(removed));

    scran::simple_pca::Options opt;
    opt.scale = scale;

    // The initial vector is slightly different when we lose a feature, so we manually force our own random initialization.
    auto raw_init = simulate_vector(nr, /* lower = */ -2, /* upper = */ 2, /* seed = */ scale + rank + 10); 
    raw_init.back() = 0;

    Eigen::VectorXd init(nr - 1);
    std::copy_n(raw_init.begin(), nr - 1, init.data());
    opt.irlba_options.initial = &init;
    auto ref = scran::simple_pca::compute(&leftovers, rank, opt);

    Eigen::VectorXd init2(nr);
    std::copy_n(raw_init.begin(), nr, init2.data());
    opt.irlba_options.initial = &init2;
    auto out = scran::simple_pca::compute(&has_zero, rank, opt);

    expect_equal_pcs(ref.components, out.components); 
    expect_equal_vectors(ref.variance_explained, out.variance_explained);
    EXPECT_FLOAT_EQ(out.total_variance, ref.total_variance);

    // Same behavior with sparse representation.
    auto sparse_zero = tatami::convert_to_compressed_sparse(&has_zero, true);
    auto spout = scran::simple_pca::compute(sparse_zero.get(), rank, opt);

    expect_equal_pcs(spout.components, out.components);
    expect_equal_vectors(spout.variance_explained, out.variance_explained);
    EXPECT_FLOAT_EQ(spout.total_variance, out.total_variance);
}

INSTANTIATE_TEST_SUITE_P(
    SimplePca,
    SimplePcaMoreTest,
    ::testing::Combine(
        ::testing::Values(false, true), // to scale or not to scale?
        ::testing::Values(2, 5, 10) // number of PCs to obtain
    )
);
