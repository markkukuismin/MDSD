# MDSD (The Mean Degree Squared Difference)

## Description

The R package `MDSD` can be used to identify hubs using solution path information. The method is tested and compatible with `glasso`, `hglasso`, and `space` R packages. See the example below.

For more details on the MDSD distance, please refer to: (link comming soon...) Supplementary data and code needed to reproduce the results reported in the article are available at: [https://github.com/markkukuismin/MDSD_supplementary](https://github.com/markkukuismin/MDSD_supplementary).

## Dependencies

Please make sure to install the following packages before using R package `MDSD`. R with version later than 4.3.1 is needed.
```r
install.packages(c("igraph", "ggplot2", "huge"))
```

## Installation

The R package `MDSD` can be installed from source files in the GitHub repository (R package `devtools` is needed):
```r
library(devtools)
install_github(repo="markkukuismin/MDSD")
```

## Main functions

* `hub_detection_mdsd`: Compute the MDSD distance from the solution path information.
* `hub_detection_plot`: Plots of the node degree with respect to the tuning parameter value or the MDSD values.
* `cor_screening`: Thresholding of the absolute values of correlation coefficients (aka. lossy screening). One can also use the function `huge` from the `huge` package [(link)](https://cran.r-project.org/web/packages/huge/index.html). I faced with some memory leak issues while using `huge` so I made this.
