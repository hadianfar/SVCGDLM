
# SVCGDLM: Spatially Varying Coefficient Generalized Distributed Lag Models 

<!-- badges: start -->
<!-- badges: end -->

SVCGDLM is an R package for fitting the Spatially Varying Coefficient Generalized Distributed Lag Model.
The package provides tools for modeling exposure–lag–response relationships that vary across space, allowing both spatial nonstationarity and spatial dependence in the estimated effects.
Parameter estimation is implemented via a Monte Carlo Expectation–Maximization (MCEM) algorithm using Polya–Gamma data augmentation, ensuring computational efficiency and scalability for large datasets.
The package includes functions for model fitting, simulation, visualization, and performance evaluation under various spatial and temporal dependency structures.


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
