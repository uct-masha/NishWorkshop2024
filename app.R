# This script creates a shiny app for a vaccine model meant to be presented
# at the NiSH Training Workshop

library(shiny)
library(bslib)
library(deSolve)
library(tidyverse)
library(sortable)
library(shinyWidgets)

# source the model in
source("model.R")

# define the ui
ui <- page_navbar(
  title = "Disease Transmission Model",
  bg = "#373A40",
  inverse = TRUE,
  theme = bs_theme(
    version = 5,
    base_font = font_google("Public Sans")
  ),
  header = tagList(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    )
  ),
  nav_spacer(),
  nav_panel(
    title = NULL,
    navset_card_pill(
      nav_panel(
        title = "Welcome",
        layout_columns(
          col_widths = c(4, 8),
          card(
            img(
              class = "logo",
              href="http://www.masha.uct.ac.za",
              src="images/logo_masha.png"
            ),
            div(
              includeMarkdown('www/markdown/about_masha.md')
            )
          ),
          accordion(
            accordion_panel(
              title = "Introduction",
              "This is an introduction to the model and the app."
            ),
            accordion_panel(
              title = "About",
              "In case we have more stuff to sa"
            ),
            accordion_panel(
              title = "How to use the App",
              "Instructions on how to use the app."
            )
          )
        )
      ),
      nav_panel(
        title = "Parameters",
        layout_column_wrap(
          class = "parameters",
          width = 1 / 3,
          heights_equal = "row",
          max_height = 400,
          gap = 100,
          numericInput(inputId = "initP", label = "Population size", value = 1000, min = 100, max = 100000, step = 100),
          sliderInput(inputId = "beta", label = "Average contact rate per year", value = 300, min = 40, max = 1000, step = 1),
          sliderInput(inputId = "nrec", label = "Average natural recovery period (days)", value = 10, min = 1, max = 20, step = 1),
          sliderInput(inputId = "lifeexp", label = "Average life expectancy (days)", value = 100, min = 80, max = 150, step = 1),
          sliderInput(inputId = "lossimm", label = "Average loss of immunity period (days)", value = 30, min = 1, max = 365, step = 1),
          sliderInput(inputId = "beta1", label = "Relative amplitude of seasonal forcing", value = 0.4, min = 0, max = 1, step = 0.01),
          sliderInput(inputId = "propv", label = "Vaccination Coverage", value = 20, min = 0, max = 100, step = 10, post = "%"),
          sliderInput(inputId = "seekvac", label = "Time to seek vaccination", value = 7, min = 1, max = 10, step = 1),
          sliderInput(inputId = "vacstart", label = "Vaccination start time", value = 10, min = 0, max = 20, step = 0.25)
        )
      ),
      nav_panel(
        title = "Scenarios",
        h4(icon("cubes"), "Build a Package of Interventions"),
        fluidRow(
          column(
            width = 6,
            source("www/R/ui/intervention_first_vaccine.R", local = TRUE)$value,
          ),
          column(
            width = 6,
            source("www/R/ui/intervention_second_vaccine.R", local = TRUE)$value
          )
        ),
        # Floating box to name and run a package of interventions
        div(
          class = "float-go-interv",
          div(
            class = "card",
            div(id = "intervention_card_header", class = "card-header", "Simulation of the Package"),
            div(
              id = "intervention_card_body", class = "card-body",
              layout_columns(
                class = "package_layout",
                span(
                  "Package Name:",
                  textInput("package_name", label = NULL, placeholder = "Vaccine 1", width = "250px")
                ),
                actionButton("go_interventions", "Run", class = "btn btn-success"),
                col_widths = c(9, 3)
              )
            )
          )
        )
      ),
      nav_panel(
        title = "Results",
        plotOutput(outputId = "model_plot")
      )
    )
  )
)

# define the server function
server <- function(input, output, session) {

  initial_conditions <- makeInitialConditions()
  contact <- makeContactMatrix(
    initialConditions = initial_conditions,
    total_contacts = total_contacts
  )

  # model output
  mod <- runModel(
    initialConditions = initial_conditions,
    parameters = makeParameters(),
    contact = contact
  )

  # plot the model output
  output$model_plot <- renderPlot({
    plotInc(mod |> filter(time>=2024))
  })
}

shinyApp(ui, server)
