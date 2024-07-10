# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c(
    "shiny",
    "bslib",
    "PBSddesolve",
    "dplyr",
    "forcats",
    "tidyr",
    "stringr",
    "tibble",
    "lubridate",
    "ggplot2",
    "docstring",
    "reactable",
    "shinycssloaders",
    "plotly"
  ),
  # format = "qs",
  controller = crew::crew_controller_local(workers = 6, seconds_idle = 60)
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
tar_source("packages.R")
tar_source("model.R")

# Replace the target list below with your own:
list(
  tar_target(
    name = pop_size,
    command = 50000000
  ),
  tar_target(
    name = initialConditions,
    command = rescalePop(unlist(makeInitialConditions()), pop_size)
  ),
  tar_target(
    name = contacts,
    command = makeContactMatrix(initialConditions = initialConditions)
  ),
  # Scenario: 0 vax1, 0 vax2 ####
  tar_target(
    name = mod_0_0,
    command = runModel(startyear=2000, endyear=2040,
                       initialConditions=initialConditions,
                       parameters=makeParameters(cov1=0, cov2=0),
                       contact = contacts)
  ),
  tar_target(
    name = plt_inc_0_0,
    command = plotModel(mod_0_0, Variables = c("Inc"), UsePlotly = F)
  ),
  tar_target(
    name = plt_prot_0_0,
    command = plotModel(mod_0_0, Variables = c("Vax"), UsePlotly = F)
  ),
  # Scenario: 85 vax1, 0 vax2 ####
  tar_target(
    name = mod_85_0,
    command = runModel(startyear=2000, endyear=2040,
                       initialConditions=initialConditions,
                       parameters=makeParameters(cov1=0.85, cov2=0),
                       contact = contacts)
  ),
  tar_target(
    name = plt_inc_85_0,
    command = plotModel(mod_85_0, Variables = c("Inc"), UsePlotly = F)
  ),
  tar_target(
    name = plt_prot_85_0,
    command = plotModel(mod_85_0, Variables = c("Vax"), UsePlotly = F)
  ),
  # Scenario: 85 vax1, 85 vax2 ####
  tar_target(
    name = mod_85_85,
    command = runModel(startyear=2000, endyear=2040,
                       initialConditions=initialConditions,
                       parameters=makeParameters(cov1=0.85, cov2=0.85),
                       contact = contacts)
  ),
  tar_target(
    name = plt_inc_85_85,
    command = plotModel(mod_85_85, Variables = c("Inc"), UsePlotly = F)
  ),
  tar_target(
    name = plt_prot_85_85,
    command = plotModel(mod_85_85, Variables = c("Vax"), UsePlotly = F)
  )
)
