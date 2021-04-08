# biokit  <img src="pics/icon.png" align="right" height="200" />
[![Build Status](https://travis-ci.com/martingarridorc/biokit.svg?branch=master)](https://travis-ci.com/martingarridorc/biokit)
[![CodeCov](https://codecov.io/gh/martingarridorc/biokit/branch/master/graph/badge.svg)](https://codecov.io/gh/martingarridorc/biokit/)
[![License](https://img.shields.io/github/license/martingarridorc/biokit.svg?color=yellow)](https://github.com/martingarridorc/biokit/blob/master/LICENSE)
![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/martingarridorc/biokit?include_prereleases)

This package is a toolkit that can be concieved as a wrapper for functions and utilities that I use repeatedly across projects and collaborations involving omics data analysis. Particularly, it makes use of the core functions from packages such as [limma](https://bioconductor.org/packages/release/bioc/html/limma.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html) to automate processes as the statistical comparison and functional analysis of omics data.

## Installation

The package can be installed through the `install_github()` function from [devtools](https://cran.r-project.org/web/packages/devtools/index.html).

```
devtools::install_github(repo = "https://github.com/martingarridorc/biokit")
```

## Acknowledgments

[Carlos Loucera](https://github.com/loucerac) and his mathematical wisdom.
