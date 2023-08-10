# Load rhandontable
library("rootmech")
library("rhandsontable")

# Initial data
data_initial <- data.frame(
  dr = seq(2, 10, l = 10),
  epsru = runif(10, 0.2, 0.3),
  tru = runif(10, 5, 8)
)
