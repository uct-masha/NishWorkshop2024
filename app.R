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
    title = "Home",
    navset_card_underline(
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
              title = "About the model",
              "Details about the model"
            ),
            accordion_panel(
              title = "How to use the App",
              "Instructions on how to use the app."
            )
          )
        )
      ),
      nav_panel(
        title = "Run Model",
        accordion(
          open = TRUE,
          fluidRow(
            column(3, fluidRow(
              column(12,
                     # sliderInput(inputId = "beta", label = "Probability of transmission", value = 10, min = 0, max = 100, post = "%"),
                     # year to start the vaccination 2024:2030
                     numericInput(inputId = "yearStart", label = "Year to start vaccination", value = 2024, min = 2024, max = 2030),
                     sliderInput(inputId = "cov1", label = "Vaccine 1 coverage", value = 65, min = 0, max = 100, post = "%"),
                     sliderInput(inputId = "cov2", label = "Vaccine 2 coverage", value = 50, min = 0, max = 100, post = "%"),
                     sliderInput(inputId = "pt", label = "Probability of seeking treatment", value = 100, min = 0, max = 100, post = "%"),
                     # numericInput(inputId = "rs", label = "Incubation rate", value = 1, min = 1, max = 10),
                     # numericInput(inputId = "rr", label = "Natural recovery rate", value = 1, min = 1, max = 10),
                     # sliderInput(inputId = "ve1", label = "Vaccine efficacy", value = 100, min = 0, max = 100, post = "%"),
                     # numericInput(inputId = "rt", label = "Treatment seeking rate", value = 1, min = 1, max = 10),
                     # numericInput(inputId = "rtr", label = "Treatment recovery rate", value = 1, min = 1, max = 10),
                     actionButton(inputId = "runModel", label = "Run Model", class = "btn btn-success")
               )
            )),
            column(12-3, fluidRow(
              card(
                full_screen = TRUE,
                plotOutput(outputId = "model_plot_incidence")
              ),
              card(
                full_screen = TRUE,
                plotOutput(outputId = "model_plot_treatment")
              )
            ))
          )
        )
      )
    )
  ),
  nav_menu(
    title = "Links",
    align = "right",
    nav_item(tags$a(
      icon("github"),
      "Code",
      href = "https://github.com/uct-masha/NishWorkshop2024",
      target = "_blank"
    )),
    nav_item(tags$a(
      icon("globe"),
      "MASHA Website",
      href = "https://science.uct.ac.za/masha",
      target = "_blank"
    )),
    nav_item(tags$a(
      icon("globe"),
      "NISH Website",
      href = "https://health.uct.ac.za/nish",
      target = "_blank"
    ))
  )
)

# global parameters
initial_conditions <- makeInitialConditions()
contact <- makeContactMatrix(
  initialConditions = initial_conditions,
  total_contacts = total_contacts
)
# define the server function
server <- function(input, output, session) {
  params <- reactive({
    makeParameters(
                   timevax = input$yearStart,
                   # beta = input$beta/100,
                   cov1 = input$cov1/100,
                   cov2 = input$cov2/100,
                   pt = input$pt/100
                   # rs = input$rs,
                   # rr = input$rr,
                   # ve1 = input$ve1/100,
                   # rt = input$rt,
                   # rtr = input$rtr
                   )
  })

  modelOutput <- reactiveVal()

  observeEvent(input$runModel, {
    # model output
    mo <- runModel(
      initialConditions = initial_conditions,
      parameters = params(),
      contact = contact
    )
    modelOutput(mo)
  })

  # plot the model output: Incidence
  output$model_plot_incidence <- renderPlot({
    plotInc(modelOutput() |> filter(time>=2024))
  })

  # plot the model output: Treatment
  output$model_plot_treatment <- renderPlot({
    plotTr(modelOutput() |> filter(time>=2024))
  })
}

shinyApp(ui, server)
