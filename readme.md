# simphony

[![CircleCI](https://circleci.com/gh/hugheylab/simphony.svg?style=shield)](https://circleci.com/gh/hugheylab/simphony)
[![codecov](https://codecov.io/gh/hugheylab/simphony/branch/master/graph/badge.svg)](https://codecov.io/gh/hugheylab/simphony)

`simphony` simulates large-scale, rhythmic data, including transcriptome data and behavioral activity data. For technical details on how we designed and validated `simphony`, check out the [paper](https://doi.org/10.7717/peerj.6985) and the [accompanying results](https://doi.org/10.6084/m9.figshare.7441355).

## Installation

If you use RStudio, go to Tools -> Global Options... -> Packages -> Add... (under Secondary repositories), then enter:

- Name: hugheylab
- Url: https://hugheylab.github.io/drat/

You only have to do this once. Then you can install or update the package by entering:

```R
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('simphony')
```

Alternatively, you can install or update the package by entering:

```R
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('simphony', site_repository = 'https://hugheylab.github.io/drat/')
```

## Usage

For an introduction to the package, read the [introductory vignette](https://simphony.hugheylab.org/articles/introduction.html). To explore more of simphony's capabilities, check out the [example vignette](https://simphony.hugheylab.org/articles/examples.html). For more details, check out the [reference documentation](https://simphony.hugheylab.org/reference/index.html).
