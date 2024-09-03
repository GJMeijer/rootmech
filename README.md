# rootmech

Power law fitting of root diameter versus root biomechanical properties.

v0.1.0 - September 2024 - Gerrit Meijer (<gjm36@bath.ac.uk>)


## Running the app online

An online app that can be used to fit user-defined data, as well as browse power law fits for a large number of plant species, is currently hosted through Shinyapps.io: <https://gjmeijer.shinyapps.io/rootmech/>


## Package and installation

This R package contains:

* an interactive app (using R Shiny), allowing you to:
  * fit custom data using various power law fitting methods in a user-friendly way
  * browse power law fit parameters for a large variety of plant species
* all R functions required to generate power law fits

To install this package on your local machine:

1. Open R. (If not installed, I recommend installing RStudio (free software) that can be downloaded from https://www.rstudio.com/products/rstudio/download/)
2. If not already installed, install the `devtools` package by typing `install.packages("devtools")` in the R console. This package allows you to interact with R packages hosted on GitHub (among many other things).
3. Install the `rootmech` package by typing `devtools::install_github("GJMeijer/rootmech")` in the R console

After installation, load the (i.e. made available for use in your current R session) by typing `library("rootmech")`.


## Running the app offline

To run the included app on your local machine, type `run_app()` in the R console after installing and loading the package. This will open the app in your default browser. This app relies on R continuing to run in the background for any computations.

Instructions for how to use the app can be found within the app itself. The app is constructed using the R package `shiny`.


## Using calculations functions

A wrapper function `powerlaw_fit()` is provided that can generate power law fits for any of the methods described in Meijer (2024). For more information about this function and its input and output arguments, use the R command `?powerlaw_fit()`.

For example, to generate a power law fit using the Weibull model a set of diameter (`x`) and tensile strength data (`y`), execute the code `powerlaw_fit(x, y, model = "weibull")`. When using custom weightings for each observation, stored in a vector `w`, execture `powerlaw_fit(x, y, model = "weibull", weights = w)` instead.

## Detailed function documentation

### Function documentation

For each of the functions in the package, a help page can be opened by executing `?function_name` in R, where `function_name` is the name of the function you want to see the documentation of. A list of all functions within 
this package can be requested by typing `library(help = "rootmech")` into the R console.

### Vignettes

Key functions are accompanied by vignettes. These include the mathematical 
derivations behind each of the fitting methods and algorithms. Execute `browseVignettes("rootmech")` 
to obtain a list with all vignettes and browse their contents (documentation in HTML). Execture `RShowDoc("name", package = "rootmech")` to open a specific vignette with the a specific name (`name`)  directly. 

Available vignettes describing power law fitting are :
* `powerlaw_gamma`: power law fitting with gamma distribution
* `powerlaw_gumbel`: power law fitting with Gumbel distribution
* `powerlaw_logistic`: power law fitting with logistic distribution
* `powerlaw_lognormal`: power law fitting with lognormal distribution
* `powerlaw_lognormal_uncorrected`: power law fitting with lognormal distribution, excluding the correction for the mean
* `powerlaw_normal`: power law fitting with normal distribution
* `powerlaw_uniform`: power law fitting with uniform distribution
* `powerlaw_weibull` power law fitting with Weibull distribution

Further vignettes are available:
* `gumbel`: Fitting a Gumbel probability distribution 
* `power`: Fitting a power law probability distribution
* `weibull`: Fitting a Weibull probability distribution


## Copyright and Licence

Copyright (C) 2024 University of Bath

This program is free software: you can redistribute it and/or modify it under the terms of the MIT license.
