# simphony
[![CircleCI](https://circleci.com/gh/hugheylab/simphony.svg?style=shield)](https://circleci.com/gh/hugheylab/simphony)

`simphony` is an R package for simulating large-scale, rhythmic data, particularly rhythmic transcriptome data. For technical details on how we designed and validated `simphony`, check out the [paper](https://doi.org/10.7717/peerj.6985) and the [accompanying results](https://doi.org/10.6084/m9.figshare.7441355).

## Installation
First install drat.
```R
install.packages('drat')
```

Then add the following line to your `.Rprofile` file (typically located at "~/.Rprofile"), which gets run every time R starts. See [here](https://csgillespie.github.io/efficientR/3-3-r-startup.html#r-startup) for details.
```R
drat::addRepo('hugheylab')
```
Alternatively, you can just enter the above command into the R console each time you want to install or update the package.

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
