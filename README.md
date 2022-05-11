# simphony

[![check-deploy](https://github.com/hugheylab/simphony/workflows/check-deploy/badge.svg)](https://github.com/hugheylab/simphony/actions)
[![codecov](https://codecov.io/gh/hugheylab/simphony/branch/master/graph/badge.svg)](https://codecov.io/gh/hugheylab/simphony)
[![Netlify Status](https://api.netlify.com/api/v1/badges/ebfe10a6-7cd8-416b-969f-eb0160f665b2/deploy-status)](https://app.netlify.com/sites/hungry-johnson-e23843/deploys)
[![CRAN Status](https://www.r-pkg.org/badges/version/simphony)](https://cran.r-project.org/package=simphony)
[![drat version](https://raw.githubusercontent.com/hugheylab/drat/gh-pages/badges/simphony_drat_badge.svg)](https://github.com/hugheylab/drat/tree/gh-pages/src/contrib)

`simphony` simulates large-scale, rhythmic data, including transcriptome data and behavioral activity data. For technical details on how we designed and validated `simphony`, check out the [paper](https://doi.org/10.7717/peerj.6985) and the [accompanying results](https://doi.org/10.6084/m9.figshare.7441355).

## Installation

### Option 1: CRAN

```r
install.packages('simphony')
```

### Option 2: Hughey Lab Drat Repository

1. Install [`BiocManager`](https://cran.r-project.org/package=BiocManager).

    ```r
    if (!requireNamespace('BiocManager', quietly = TRUE))
      install.packages('BiocManager')
    ```

1. If you use RStudio, go to Tools → Global Options... → Packages → Add... (under Secondary repositories), then enter:

    - Name: hugheylab
    - Url: https://hugheylab.github.io/drat/

    You only have to do this once. Then you can install or update the package by entering:

    ```r
    BiocManager::install('simphony')
    ```

    Alternatively, you can install or update the package by entering:

    ```r
    BiocManager::install('simphony', site_repository = 'https://hugheylab.github.io/drat/')
    ```

## Usage

For an introduction to the package, read the [introductory vignette](https://simphony.hugheylab.org/articles/introduction.html). To explore more of simphony's capabilities, check out the [example vignette](https://simphony.hugheylab.org/articles/examples.html). For more details, check out the [reference documentation](https://simphony.hugheylab.org/reference/index.html).
