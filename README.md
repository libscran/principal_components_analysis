# Principal components analysis, duh

![Unit tests](https://github.com/libscran/principal_component_analysis/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/libscran/principal_component_analysis/actions/workflows/doxygenate.yaml/badge.svg)
[![Codecov](https://codecov.io/gh/libscran/principal_component_analysis/graph/badge.svg?token=qklLZtJSE9)](https://codecov.io/gh/libscran/principal_component_analysis)

## Overview

As the name suggests, this repository implements functions to perform a PCA on the gene-by-cell expression matrix,
returning low-dimensional coordinates for each cell that can be used for efficient downstream analyses, e.g., clustering, visualization.
The code itself was originally derived from the [**scran**](https://bioconductor.org/packages/scran) and [**batchelor**](https://bioconductor.org/packages/batchelor) R packages
factored out into a separate C++ library for easier re-use.

## Quick start

Given a [`tatami::Matrix`](https://github.com/tatami-inc/tatami), the `simple_pca::compute()` function will compute the PCA to obtain a low-dimensional representation of the cells:

```cpp
#include "scran/simple_pca.hpp"

tatami::Matrix<double, int>* ptr = some_data_source();

// Take the top 20 PCs:
scran::simple_pca::Options opt;
auto res = scran::simple_pca::compute(ptr, 20, opt);

res.components; // rows are PCs, columns are cells.
res.rotation; // rows are genes, columns correspond to PCs.
res.variance_explained; // one per PC, in decreasing order.
res.total_variance; // total variance in the dataset.
```

Advanced users can fiddle with the options:

```cpp
opt.scale = true;
opt.num_threads = 4;
opt.realize_matrix = false;
auto res2 = scran::simple_pca::compute(ptr, 20, opt);
```

In the presence of multiple blocks, we can perform the PCA on the residuals after regressing out the blocking factor.
This ensures that the inter-block differences do not contribute to the first few PCs, instead favoring the representation of intra-block variation.

```cpp
std::vector<int> blocks = some_blocks();

scran::blocked_pca::Options bopt;
auto bres = scran::blocked_pca::compute(ptr, blocks.data(), 20, bopt);

bres.components; // rows are PCs, columns are cells.
bres.center; // rows are blocks, columns are genes.
```

The components derived from the residuals will only be free of inter-block differences under certain conditions (equal population composition with a consistent shift between blocks).
If this is not the case, more sophisticated batch correction methods are required.
If those methods accept a low-dimensional representation for the cells as input, 
we can use `blocked_pca::compute()` to obtain an appropriate matrix that focuses on intra-block variation without making assumptions about the inter-block differences:

```cpp
bopt.components_from_residuals = false;
auto bres2 = scran::blocked_pca::compute(ptr, blocks.data(), 20, bopt);
```

Check out the [reference documentation](https://libscran.github.io/principal_component_analysis) for more details.

## Building projects

This repository is part of the broader [**libscran**](https://github.com/libscran/libscran) library,
so users are recommended to use the latter in their projects.
**libscran** developers should just use CMake with `FetchContent`:

```cmake
include(FetchContent)

FetchContent_Declare(
  scran_principal_component_analysis 
  GIT_REPOSITORY https://github.com/libscran/principal_component_analysis
  GIT_TAG master # or any version of interest
)

FetchContent_MakeAvailable(scran_principal_component_analysis)

# For executables:
target_link_libraries(myexe scran_principal_component_analysis)

# For libaries
target_link_libraries(mylib INTERFACE scran_principal_component_analysis)
```
