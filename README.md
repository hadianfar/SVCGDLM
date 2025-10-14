
# SVCGDLM: Spatially Varying Coefficient Generalized Distributed Lag Models 

<!-- badges: start -->
<!-- badges: end -->

This package implements a Monte Carlo Expectation-Maximization (MCEM) algorithm to estimate Spatially varying coefficient. The method relies using Polya-Gamma data augmentation to facilitate closed-form conditional distributions in the expectation step via Gibbs sampling. The maximization step then optimizes the expected complete-data log likelihood over both fixed and spatial parameters.



## Installation

```r
# Install the development version from GitHub
# install.packages("devtools")
devtools::install_github("hadianfar/SVCGDLM")
```

## Usage

```r
library(SVCGDLM)
# Example usage:
# result <- MCEM(mcmc_samples = 1000, y = ..., x = ..., ...)
```

## Documentation

- See [function reference](./man/) and vignettes (if available).

## Author

Ali Hadianfar

## License

MIT
