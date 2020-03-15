# MGR biokit  <img src="icon/icon.png" align="right" height="200" />
[![Build Status](https://travis-ci.com/martingarrido/biokit.svg?branch=master)](https://travis-ci.com/martingarrido/biokit)
[![CodeCov](https://codecov.io/gh/martingarrido/biokit/branch/master/graph/badge.svg)](https://codecov.io/gh/martingarrido/biokit/)

This package is my personal wrapper for several functions and utilities from different packages which I use repeatedly across projects and collaborations involving omic data analysis. Particularly, it makes use of the core functions from the packages [limma](https://bioconductor.org/packages/release/bioc/html/limma.html),  [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) to automate processes such as the statistical comparison and functional analysis of omic data wrapped in the [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) class.

## Installation

The package can be installed by running the `install_github()` function from the [devtools package](https://cran.r-project.org/web/packages/devtools/index.html).

`devtools::install_github(repo = "https://github.com/martingarrido/biokit")`
