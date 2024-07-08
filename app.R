# This script creates a shiny app for a vaccine model meant to be presented
# at the NiSH Training Workshop

#notes
# removed vaccination from R1
# introduced waning functionality in Age 3 only - not used
# adjusted pop mix to SA pop
# created population protected plot
# added a fourth age group (2-5 years)
## need to add costs to UI
## need to add a cost function to the model
## need to add table to the model output
## need to engage with costing and question structure and create rubric.

library(shiny)
library(bslib)
library(PBSddesolve)
library(dplyr)
library(forcats)
library(tidyr)
library(stringr)
library(tibble)
library(lubridate)
library(ggplot2)
library(docstring) # Lets you use ?func for functions in this file

# source the model in
source("model.R")

img_src <- as.character(paste0('data:image/png;base64,', readLines("logo_masha.txt")[[1]], collapse=""))

# define the ui ####
ui <- page_navbar(
  title = "Disease Transmission Model",
  bg = "#373A40",
  inverse = TRUE,
  theme = bs_theme(
    version = 5,
    base_font = "Roboto",
    base_font_google = NULL,
    code_font = "Roboto Mono",
    code_font_google = NULL
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
            img(src = img_src, class = "logo"),
            tags$p("The ",
                   a("Modelling and Simulation Hub, Africa (MASHA)", href="http://www.masha.uct.ac.za"),
                   " is a research group at the University of Cape Town. MASHA’s research focus is the development and application of mathematical modelling and computer simulation to predict the dynamics and control of infectious diseases to evaluate the impact of policies aimed at reducing morbidity and mortality. Based in the Faculty of Science, MASHA’s research is closely integrated with other disciplines resulting in policy-driven and impactful scientific research."
            )
          ),
          div(
            tags$h1("Disease Transmission Model"),
            tags$p("This application has been developed to exemplify how mathematical modelling can be used to provide scientific evidence to support decisions on vaccine introduction. In this application, a transmission model of a disease X has been developed for a population where the disease is currently transmitting.
            Disease X has been circulating in the population for many years. It has the following properties:"),
            tags$ul(
            tags$li("Having the disease confers life-long immunity (similar to measles)"),
            tags$li("In the absence of vaccination, there are high levels of immunity in the older populations, but the very young are left without protection."),
            tags$li("A new vaccine has been developed that is available globally."),
            tags$li("This vaccine can be delivered in 1 or 2 dose format, where the first dose has protection efficacy of 85% and having both doses confers full (100%) protection."),
            tags$li("The recommended schedule from the Global Health Authority is that the two doses are delivered at 1 year and 2 years of age respectively."),
            tags$li("The planned introduction of the vaccine is 2025."),
            tags$li("The vaccine is not yet included in the national immunisation schedule.")
                     ),
            tags$p("This app was created for the NiSH Training Workshop, which is a training workshop for health economists and modellers in South Africa. The workshop is organised by the ",
                   a("National Institute for Health Innovation (NiSH)", href="https://www.nish.ac.za"),
                   " and the ",
                   a("Modelling and Simulation Hub, Africa (MASHA)", href="http://www.masha.uct.ac.za"),
                   "."
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
                     # year to start the vaccination 2024:2040
                     p("Vaccination"),
                     numericInput(inputId = "yearStart", label = "Year to start vaccination", value = 2025, min = 2025, max = 2025),
                     sliderInput(inputId = "cov1", label = "Vaccine 1 coverage", value = 0, min = 0, max = 100, post = "%"),
                     sliderInput(inputId = "cov2", label = "Vaccine 2 coverage", value = 0, min = 0, max = 100, post = "%"),
                     sliderInput(inputId = "pt", label = "Probability of seeking treatment", value = 100, min = 0, max = 100, post = "%"),
                     p("Costs in USD"),
                     numericInput(inputId = "cvacc", label = "Cost per vaccine", value = 1.35, min = 1.35, max = 1.35),
                     numericInput(inputId = "cdel", label = "Cost per vaccine delivered", value = 3, min = 1, max = 1000),
                     numericInput(inputId = "ctrt", label = "Cost per case treated", value = 1.2, min = 1, max = 1000),
                     numericInput(inputId = "cintro", label = "Introduction cost (once-off)", value = 850,000, min = 1, max = 10000000),
                     # numericInput(inputId = "rs", label = "Incubation rate", value = 1, min = 1, max = 10),

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
                plotOutput(outputId = "model_plot_protected")
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

# global parameters ####
initial_conditions <- makeInitialConditions()
total_contacts <- matrix(c(8,  5,  3,  15,
                           5,  14, 5,  12,
                           3,  5,  20, 12,
                           15, 12, 12, 30), nrow=4)*365
contact <- makeContactMatrix(
  initialConditions = initial_conditions,
  total_contacts = total_contacts
)
mo_initial <- runModel(
  initialConditions = initial_conditions,
  parameters = makeParameters(),
  contact = contact
)

# define the server function ####
server <- function(input, output, session) {
  params <- reactive({
    makeParameters(
                   timevax = input$yearStart,
                   # beta = input$beta/100,
                   cov1 = input$cov1/100,
                   cov2 = input$cov2/100,
                   pt = input$pt/100,
                   cvacc = input$cvacc,
                   cdel = input$cdel,
                   ctrt = input$ctrt,
                   cintro = input$cintro
                   # rs = input$rs,
                   # rr = input$rr,
                   # ve1 = input$ve1/100,
                   # rt = input$rt,
                   # rtr = input$rtr
                   )
  })

  modelOutput <- reactiveVal(mo_initial)

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
    plotInc(modelOutput() |> filter(time>=2022))
  })

  # plot the model output: Treatment
  output$model_plot_protected <- renderPlot({
    plotProt(modelOutput() |> filter(time>=2022))
  })}

shinyApp(ui, server)
