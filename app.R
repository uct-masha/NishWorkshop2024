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
library(reactable)
library(shinycssloaders)

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
    ),
    shinyjs::useShinyjs(), # Required to enable shinyjs functionality
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
                     numericInput(inputId = "cintro", label = "Introduction cost (once-off)", value = 850000, min = 1, max = 10000000),
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
                card_header("Outputs are shown as the total for the full model timeframe (2025 - 2040)"),
                reactableOutput(outputId = "model_table") %>%
                  withSpinner()
              ),
              navset_card_tab(
                full_screen = TRUE,
                nav_panel(
                  title = "Incidence",
                  plotOutput(outputId = "model_plot_incidence") %>%
                    withSpinner()
                ),
                nav_panel(
                  title = "Protected",
                  plotOutput(outputId = "model_plot_protected") %>%
                    withSpinner()
                ),
                nav_panel(
                  title = "Treated",
                  plotOutput("model_plot_treated") %>%
                    withSpinner()
                )
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
    showNotification(span(h4(icon("hourglass-half"), "Model Running..."), "typically runs in 2 to 10 secs."),
      duration = NULL, type = "message", id = "model_run"
    )
    shinyjs::disable("runModel")
    # model output
    mo <- runModel(
      initialConditions = initial_conditions,
      parameters = params(),
      contact = contact
    )
    modelOutput(mo)
    removeNotification(id = "model_run", session = session)
    shinyjs::enable("runModel")
  })

  # table
  output$model_table <- renderReactable({

    tbIncAges <- modelOutput() |>
      mutate(age=str_extract(compartment,pattern="\\d+$") |> as.numeric(),
             Year=year(date_decimal(time)),
             Month=month(date_decimal(time))) |>
      filter(compartment |> str_starts("CInc")) |>
      summarise(Incidence=(last(population) - first(population)),
                .by=c(Year, age))

    tbProtAges <- modelOutput() |>
      mutate(age=str_extract(compartment,pattern="\\d+$") |> as.numeric(),
             Year=year(date_decimal(time)),
             Month=month(date_decimal(time))) |>
      filter(compartment %in% c("R1", "R2", "R3", "R4",  "VA2", "VA3", "VB3", "VA4", "VB4"),
             Year==(time)) |>
      summarise(Protected=sum(population),
                .by=c(Year, age))

    tbTrAges <- modelOutput() |>
      mutate(age=str_extract(compartment,pattern="\\d+$") |> as.numeric(),
             Year=year(date_decimal(time)),
             Month=month(date_decimal(time))) |>
      filter(compartment |> str_starts("CTr")) |>
      summarise(Treatment=last(population) - first(population),
                .by=c(Year, age))


    result <- tbIncAges %>%
      left_join(tbProtAges, by = join_by(Year, age)) %>%
      left_join(tbTrAges, by = join_by(Year, age)) %>%
      filter(between(Year, 2025, 2040)) %>%
      summarise(
        total_incidence = sum(Incidence),
        total_protected = sum(Protected),
        total_treated = sum(Treatment),
        .by = age
      )
    reactable(
      result,
      defaultColDef = colDef(
        format = colFormat(digits = 0, separators = TRUE)
      ),
      columns = list(
        total_incidence = colDef(
          name = "Total Incidence"
        ),
        total_protected = colDef(
          name = "Total Protected"
        ),
        total_treated = colDef(
          name = "Total Treated"
        )
      )
    )
  })

  # plot the model output: Incidence
  output$model_plot_incidence <- renderPlot({
    plotInc(modelOutput() |> filter(time>=2022))
  })

  # plot the model output: Protected
  output$model_plot_protected <- renderPlot({
    plotProt(modelOutput() |> filter(time>=2022))
  })

  # plot the model output: Treatment
  output$model_plot_treated <- renderPlot({
    plotTr(modelOutput() |> filter(time>=2022))
  })
}

  # plot the model output: Treatment

shinyApp(ui, server)
