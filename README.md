# reflbrown

## Overview

R package reflbrown provides functions for

-   random number generation of first hitting time distribution of the reflected Brownian motion;

-   simulating data including recurrent event time and covariates

-   Fitting Independent-frailty model


## Installation

### Install dependences:

install.packages(c(
  "Rcpp",
  "copula",
  "ggplot2",
  "reda",
  "Formula",
  "loo",
  "nimble"
))


### Installation


User can install from the Github

```r
install.packages("remotes")
remotes::install_github("YingfaX/reflbrown")
```

Or download the package and install with 
```r 
install.packages("reflbrown_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Citation

**BibTeX Entry:**

``` bibtex
@misc{reflbrown-package,
    title = {{reflbrown}: {R}ecurrent Event Modeling based on the First Hitting Time
    of {R}eflected Brownian Motion},
    author = {Yingfa Xie and Jun Yan},
    year = {2024},
    howpublished = {https://github.com/YingfaX/reflbrown},
    note = {{R} package version 0.1.0}
  }
```

``` bibtex
@article{xie2025recurrent,
  title={Recurrent events modeling based on a reflected Brownian motion with application to hypoglycemia},
  author={Xie, Yingfa and Fu, Haoda and Huang, Yuan and Pozdnyakov, Vladimir and Yan, Jun},
  journal={Biostatistics},
  volume={26},
  number={1},
  pages={kxae053},
  year={2025},
  publisher={Oxford University Press}
}
```
