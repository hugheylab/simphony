# simphony
[![CircleCI](https://circleci.com/gh/hugheylab/simphony.svg?style=shield)](https://circleci.com/gh/hugheylab/simphony)

`simphony` is an R package for simulating large-scale, rhythmic data, including transcriptome data and behavioral activity data. For technical details on how we designed and validated `simphony`, check out the [paper](https://doi.org/10.7717/peerj.6985) and the [accompanying results](https://doi.org/10.6084/m9.figshare.7441355).

## Installation

First add the hugheylab repository to your repos. There are multiple ways to do this.

If you use RStudio, go to Tools -> Global Options... -> Packages -> Add... (under Secondary repositories), then enter the following values.

- Name: hugheylab
- Url: https://hugheylab.github.io/drat/

You only have to do this once.

Alternatively, you can enter the following command each time you want to install or update the package.
```R
options(repos = c(getOption('repos'), 'https://hugheylab.github.io/drat/'))
```

Now you can install the package.
```R
install.packages('simphony', type = 'source')
```
You can update the package using `update.packages()`.

## Docker
You can also use a pre-built [docker image](https://hub.docker.com/r/hugheylab/hugheyverse), which has all dependencies installed.
```bash
docker pull hugheylab/hugheyverse
```

## Getting started
Check out the vignette and the documentation.
```R
browseVignettes('simphony')
?simphony::simphony
```
