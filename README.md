# cheapknockoff

The `cheapknockoff` package contains implementation of a multiple-knockoff based procedure for performing cost-conscious feature selection when variables have costs.
The details of the method can be found in 
[Yu, Witten, and Bien (2019) *Controlling Costs: Feature Selection on a Budget*](https://arxiv.org/abs/1910.03627).

To install `cheapknockoff` from [github](http://github.com), type in R console
```R
devtools::install_github("hugogogo/cheapknockoff", build_vignettes = TRUE)
```
Note that the installation above requires using R package [devtools](https://CRAN.R-project.org/package=devtools)
(which can be installed using `install.packages("devtools")`).

Please check the accompanying vignette on how to use the `cheapknockoff` package. To read vignette, after installing the package, type in R console
```R
browseVignettes("cheapknockoff")
```
