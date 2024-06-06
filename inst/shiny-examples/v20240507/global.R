## INITIATE ####
# load packages
library("tidyverse")
library("plotly")
library("rhandsontable")

## SETTINGS ####
# settings
digits <- 4  # number of significant digits to use for fit parameters displayed
# names
user_data_name <- "<User data - copy-paste into table below>"
# default custom data
custom_default <- data.frame(dr = as.double(rep(NA, 3)), tru = as.double(rep(NA, 3)))

## LOAD DATA ####
# root
root <- "./inst/shiny-examples/v20240507/data/"
root <- "./data/"
# load files
dm <- readr::read_csv(paste0(root, "measurements.csv")) %>%
  rename(dr = `diameter [mm]`, tru = `tensile strength [MPa]`)
ds <- readr::read_csv(paste0(root, "species.csv"))
dd <- readr::read_csv(paste0(root, "datasets.csv"))
dl <- readr::read_csv(paste0(root, "sources.csv"))
# generate dataset properties - all in single dataframe
df <- dd %>%
  select(dataset_id, source_id, species_id, notes) %>%
  left_join(ds %>% select(species_id, species), by = "species_id") %>%
  left_join(dl, by = "source_id") %>%
  mutate(
    notes_text = ifelse(
      is.na(notes),
      "",
      paste0(" (", notes, ")")
    ),
    label = paste0(species, notes_text, "; ", reference),
    link_label = ifelse(
      is.na(doi),
      url,
      paste0("DOI: ", doi)
    ),
    link = ifelse(
      is.na(doi),
      url,
      paste0("https://www.doi.org/", doi)
    )
  )
# dataset options
dataset_options <- c(user_data_name, sort(unique(df$label)))

## FIT OPTIONS ####
# fit options
df_fittype <- data.frame(
  name = c(
    "Gamma",
    "Gumbel",
    "Logistic",
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
    "lognormal",
    "normal_strength",
    "normal_force",
    "normal_scaled",
    "uniform",
    "weibull"
  )
)
