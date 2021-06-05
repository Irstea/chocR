chocR: a R package implementing choc analysis
==================================================

<img src="man/figures/logo.png" align="right" width="220" />

---
License: BSD3

# Installation #
The easiest solution is to use the 'devtools' packages, this will installed the package and all required dependencies. To use 'devtools' , you need to install first:
* on Windows: [Rtools](http://cran.r-project.org/bin/windows/Rtools/)  
* on Mac: [Xcode command line tools](https://developer.apple.com/downloads)  
* on Linux: the R development package, usually called r-devel or r-base-dev  
  
On Mac, you may also need to add specific tools for Rcpp and RcppArmadillo, have a look [there](https://thecoatlessprofessor.com/programming/cpp/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/) to see the solution depending on your MAC OS version
To generate the vignette, you also need to have pandoc and pandoc-citeproc installed. Instructions can be found [here](https://pandoc.org/installing.html).    
  
Then, on a R console:

    > install.packages("devtools")
    > library(devtools)
    > install_github("Irstea/chocR",build=TRUE,build_opts = c("--no-resave-data","--no-manual"),build_vignettes = TRUE, dependencies = TRUE)

# Usage #
A vignette is included in the package to explain the usage.  

# Bug Reporting #
Please, report bug on the [github site](https://github.com/Irstea/chocR/issues).
