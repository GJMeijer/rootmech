# rootmech

Power law fitting of root diameter versus root biomechanical properties.

v0.1.0 - June 2024 - Gerrit Meijer (<gjm36@bath.ac.uk>)


## Running the app online

An online app that can be used to fit user-defined data is currently hosted through Shinyapps.io: <https://gjmeijer.shinyapps.io/rootmech/>


## Package and installation

This R package contains:

* an interactive app (using R Shiny), allowing you to fit custom data using 
various power law fitting methods in a user-friendly way
* all functions required to generate fits

To install this package on your local machine:

1. Open R. (If not installed, I recommend installing RStudio (free software) that can be downloaded from https://www.rstudio.com/products/rstudio/download/)
2. If not already installed, install the `devtools` package by typing `install.packages("devtools")` in the R console. This package allows you to interact with R packages hosted on GitHub (among many other things)
3. Install the `rootmech` package by typing `devtools::install_github("GJMeijer/rootmech")` in the R console


## Running the app offline

To run the included app on your local machine, type `rootmech::run_app()` in the R console after installation of the package. This will open the app in your default browser. This app relies on R continuing to run in the background for any computations.

Instructions for how to use the app can be found within the app itself. The app is constructed using the R package `shiny`.


## Using calculations functions

A wrapper function `powerlaw_fit()` is provided that can generate power law fits for any of the methods described in Meijer (2024). For more information about this function and its input and output arguments, type `?powerlaw_fit()` (or `?rootmech::powerlaw_fit()` if the `rootmech` packages is not currently loaded) into your R console.


## Detailed function documentation

### Function documentation

For each of the functions in the package, a help page can be opened by typing `?function_name` in the R console (or `?rootmech::function_name`, in case the `rootmech` package is has not currently loaded), where `function_name` is the name of the function you want to see the documentation of. A list of all functions within 
this package can be requested by typing `library(help = "rootmech")` into the R console.

### Vignettes

Key functions are accompagnied by vignettes which include the mathematical 
derivation behind each of the fitting methods. Type `browseVignettes("rootmech")` 
in the console to obtain a list with all vignettes and browse their contents (documentation in HTML). Type `RShowDoc("name", package = "rootmech")` to open a specific vignette with the name `"name"` directly. Vignette names correspond with function names.

## Copyright and Licence

Copyright (C) 2024 University of Bath

This program is free software: you can redistribute it and/or modify it under the terms of the MIT license.
