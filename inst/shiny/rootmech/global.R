## INITIATE ####

## SETTINGS ####
# settings
digits <- 4  # number of significant digits to use for fit parameters displayed

## DATA ####
# root data - use Zea Mays data from Meijer et al. 2024
df_maize <- data.frame(
  dr = c(
    0.65, 1.02, 1.17, 1.19, 1.37, 1.21, 2.03, 1.95, 1.81, 1.90, 3.80, 1.41,
    1.55, 2.45, 1.43, 1.03, 1.66, 1.00, 0.81, 1.19, 1.25, 1.20, 0.97, 1.94,
    1.42, 1.96, 1.31, 1.31, 2.74, 2.03, 1.94, 2.81, 2.31, 3.89, 2.59, 1.02,
    1.83, 2.02, 1.17, 1.23, 1.03, 1.53, 1.33, 1.92, 1.94, 1.99
  ),
  tru = c(
    10.918217, 9.977626,  7.504198,  12.294525, 8.519017,  7.210183,  14.574791,
    8.458128,  7.338782,  9.224462,  3.016449,  9.878638,  7.510126,  5.603739,
    7.666584,  5.207453,  3.386403,  9.632057,  16.700959, 5.930576,  8.574097,
    13.898647, 17.870551, 6.323233,  13.211660, 11.053679, 19.390546, 11.507456,
    9.166015,  7.967750,  12.667103, 9.817011,  8.787479,  3.725034,  4.621593,
    19.988294, 4.405335,  3.292624,  6.737780,  10.280848, 13.375676, 7.967198,
    9.031227,  10.549196, 6.272149,  11.754341
  ),
  weight = 1
)
df_maize <- df_maize[order(df_maize$dr), ]

## FIT OPTIONS ####
# fit options
df_fittype <- data.frame(
  name = c(
    "Gamma",
    "Gumbel",
    "Logistic",
    "Lognormal (uncorrected)",
    "Lognormal (corrected)",
    "Normal (strength)",
    "Normal (force)",
    "Normal (scaled)",
    "Uniform",
    "Weibull"
  ),
  label = c(
    "gamma",
    "gumbel",
    "logistic",
    "lognormal_uncorrected",
    "lognormal",
    "normal_strength",
    "normal_force",
    "normal_scaled",
    "uniform",
    "weibull"
  )
)

# Load fit data
if (FALSE) {
  df_fits <- readr::read_csv("./www/20240509_results_bestfit.csv")
  # merge species
  df_spec <- readr::read_csv("./www/species.csv")
  df_spec <- df_spec[, c("species", "common_name", "family")]
  df_fits <- dplyr::left_join(df_fits, df_spec, by = "species")
  readr::write_csv(df_fits, "./www/powerlawfits.csv")
}
df_fits <- readr::read_csv("./www/powerlawfits.csv")
df_fits <- df_fits[order(df_fits$species, df_fits$reference, df_fits$notes), ]
# add lognormal (uncorrected data)
if (FALSE) {
  df_fits_add <- df_fits[df_fits$model == "lognormal", ]
  df_fits_add$`t_ru0 [MPa]` <- df_fits_add$`t_ru0 [MPa]`/exp(0.5*df_fits_add$sdlog^2)
  df_fits_add$model <- "lognormal_uncorrected"
  df_fits <- dplyr::bind_rows(df_fits, df_fits_add)
  df_fits$model[df_fits$model=="lognormal_corrected"] <- "lognormal"
  readr::write_csv(df_fits, "bla.csv")
}

df_fits$label <- ifelse(
  is.na(df_fits$notes),
  paste0(df_fits$species, ", ", df_fits$reference),
  paste0(df_fits$species, ", ", df_fits$reference, ", ", df_fits$notes)
)
fit_opts <- unique(df_fits[, c("functional_group2", "family", "species")])
fit_opts <- fit_opts[order(fit_opts$functional_group2, fit_opts$family, fit_opts$species), ]

