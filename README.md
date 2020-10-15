
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multibridge: Evaluating Multinomial Order Restrictions with Bridge Sampling

<!-- badges: start -->

[![CRAN/METACRAN](https://img.shields.io/cran/v/multibridge?label=CRAN&logo=r)](https://cran.r-project.org/web/packages/multibridge/index.html)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![GitHub last commit
(devel)](https://img.shields.io/github/last-commit/crsh/multibridge/devel?label=Last%20commit&logo=github&logoColor=%23FFF)
![Travis build
status](https://img.shields.io/travis/crsh/multibridge?label=Build&logo=travis-ci&logoColor=%23FFF)
[![codecov](https://codecov.io/gh/crsh/multibridge/branch/master/graph/badge.svg)](https://codecov.io/gh/crsh/multibridge)
[![GitHub bug
issues](https://img.shields.io/github/issues/crsh/multibridge/bug?label=Bugs&logo=github&logoColor=%23FFF)](https://github.com/crsh/multibridge/issues?q=is%3Aopen+is%3Aissue+label%3Abug)
<!-- badges: end -->

Evaluates hypotheses concerning the distribution of multinomial
proportions using bridge sampling. The bridge sampling routine is able
to compute Bayes factors for hypotheses that entail inequality
constraints, equality constraints, free parameters, and a mix of all
three. These hypotheses are tested against the encompassing hypothesis,
that all parameters vary freely.

## Installation

You can install the released version of `multibridge` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("multibridge")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("ASarafoglou/multibridge")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("multibridge")
## basic example code
```

## Package dependencies

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />
