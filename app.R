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
library(plotly)
library(htmltools)

# source the model in
source("model.R")

img_src <- as.character(paste0('data:image/png;base64,', readLines("logo_masha.txt")[[1]], collapse=""))

# define the ui ####
ui <- page_navbar(
  title = "VaxSim",
  navbar_options = navbar_options(
    bg = "#373A40",
    theme = "auto"
  ),
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
            tags$p("This app was originally created for the 3rd Annual Vaccinology Course for NITAGs (AVCN) 2024'; a training workshop for National Immunisation Technical Advisory Groups in Africa. The workshop is organised by the ",
                   a("Vaccines for Africa Initiative Nitag Support Hub ", href="https://health.uct.ac.za/vacfa/nish"),
                   " with modelling support from the ",
                   a("Modelling and Simulation Hub, Africa (MASHA)", href="http://www.masha.uct.ac.za"),
                   "."
            ),
            tags$p(),
            tags$p("The ",
                   a("Modelling and Simulation Hub, Africa (MASHA)", href="http://www.masha.uct.ac.za"),
                   " is a research group at the University of Cape Town. MASHA’s research focus is the development and application of mathematical modelling and computer simulation to predict the dynamics and control of infectious diseases to evaluate the impact of policies aimed at reducing morbidity and mortality. Based in the Faculty of Science, MASHA’s research is closely integrated with other disciplines resulting in policy-driven and impactful scientific research."),
          ),
          div(
            tags$h3("About VaxSim"),
            tags$p("VaxSim is computer simulation that has been developed to exemplify how mathematical modelling can be used to provide scientific evidence to support decisions on vaccine introduction. In this simulation, a transmission model of a Disease X has been developed for a population where the disease is currently transmitting."),
            tags$h3("Population and disease characteristics"),
            tags$p("The study population is a country of 50 million people in 2024. Disease X has been circulating in the population for many years. The disease has the following properties:"),
            tags$ul(
            tags$li("The Reproductive Number (R0) of the disease is 12, with a latency period of 10 days and infectious period of 30 days"),
            tags$li("Having the disease confers life-long immunity. So, in the absence of vaccination, there are high levels of immunity in the older populations, but the very young are left without protection."),
            tags$li("All infections are symptomatic. Diagnosis and treatment can reduce the infectious period to 10 days. For example, if 70% of cases are diagnosed and treated, the effective Reproductive Number (Rₑ) is reduced to 6."),
            ),
            tags$p("A new vaccine has been developed that is available globally."),
            tags$ul(
            tags$li("This vaccine can be delivered in 1 or 2 dose format, where the first dose has protection efficacy of 85% and having both doses confers full (100%) protection."),
            tags$li("The recommended schedule from the Global Health Authority is that the two doses are delivered at 1 year and 2 years of age respectively."),
            tags$li("The planned introduction of the vaccine is 2025."),
            tags$li("The vaccine is not yet included in the national immunisation schedule.")
                     ),
            tags$h3("How to use the App"),
            tags$p("Use the sliders and input boxes on the left panel to select intended operational coverage of 1 or both doses of vaccine, and the costs of introducing the vaccine. Once you have made your selection, click the run button to run the transmission model to see the impact of vaccination on disease incidence, cases treated and the protection levels in the population. Use the summary table provided to assess the total cases and costs of vaccination for every scenario you create. You can use this output to generate the cases averted by vaccination and make a ratio of cost per case averted.")
          )
        )
      ),
      nav_panel(
        title = "Run Model",
          fluidRow(
            column(3, fluidRow(
            column(
              width = 12,
              accordion(
                open = 1,
                accordion_panel(
                  title = "Vaccination",
                  numericInput(inputId = "yearStart", label = "Year to start vaccination", value = 2025, min = 2025, max = 2025),
                  sliderInput(inputId = "cov1", label = "Vaccine 1 coverage", value = 0, min = 0, max = 100, post = "%"),
                  sliderInput(inputId = "cov2", label = "Vaccine 2 coverage", value = 0, min = 0, max = 100, post = "%"),
                  sliderInput(inputId = "pt", label = "Probability of seeking treatment", value = 70, min = 0, max = 100, post = "%")
                ),
                accordion_panel(
                  title = "Costs in USD",
                  numericInput(inputId = "cvacc", label = "Cost per vaccine", value = 1.35, min = 0, max = 1.35),
                  numericInput(inputId = "cdel", label = "Cost per vaccine delivered", value = 0.2, min = 0, max = 1000),
                  numericInput(inputId = "ctrt", label = "Cost per case treated", value = 0.75, min = 1, max = 1000),
                  numericInput(inputId = "cintro", label = "Introduction cost (once-off)", value = 500000, min = 1, max = 10000000),
                )
              ),
              br(),
              layout_columns(actionButton(inputId = "runModel", label = "Run Model", class = "btn-success"))
              )


            )),
            column(12-3, fluidRow(
              card(
                card_header("Outputs are shown as the total for the full model timeframe (2025 - 2040). All costs are discounted at 3% p.a."),
                reactableOutput(outputId = "model_table") %>%
                  withSpinner()
              ),
              navset_card_tab(
                full_screen = TRUE,
                nav_panel(
                  title = "Incidence",
                  plotlyOutput(outputId = "model_plot_incidence") %>%
                    withSpinner()
                ),
                nav_panel(
                  title = "Protected",
                  plotlyOutput(outputId = "model_plot_protected") %>%
                    withSpinner()
                ),
                nav_panel(
                  title = "Treated",
                  plotlyOutput("model_plot_treated") %>%
                    withSpinner()
                )
              )
            ))
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
      icon("chart-line"),
      "What is modelling?",
      href = "www.masha.uct.ac.za/resources/what-model",
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
initial_conditions <- rescalePop(unlist(makeInitialConditions()), 50000000)
contact <- makeContactMatrix(initialConditions = initial_conditions)
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
  cost_table <- reactiveVal(costing(mod = mo_initial, params = makeParameters()))

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
    costs <- costing(mod = mo, params = params())
    modelOutput(mo)
    cost_table(costs)
    removeNotification(id = "model_run", session = session)
    shinyjs::enable("runModel")
  })

  # cost table
  output$model_table <- renderReactable({

    cost_tbl <- cost_table() %>%
      pivot_wider(
        names_from = "variable",
        values_from = "Value"
      ) %>%
      mutate(
        discounting_factor = (1 + 0.03) ^ (2025 - Year),
        total_costs_discounted = CostTot * discounting_factor
      ) %>%
      select(!discounting_factor) %>%
      summarise(
        total_incidence = sum(Inc),
        total_costs = sum(total_costs_discounted),
        .by = age
      )

    costs_total <- sum(cost_tbl$total_costs) + as.numeric(input$cintro)
    total_cost_per_case = costs_total / sum(cost_tbl$total_incidence)

    reactable(
      cost_tbl %>% mutate(
        age = factor(age),
        age = fct_recode(age, "0-1yr" = "1", "1-2yrs" = "2", "2-5yrs" = "3", ">5yrs" = "4"),
        cost_per_case = total_costs / total_incidence
      ),
      defaultColDef = colDef(
        format = colFormat(digits = 0, separators = TRUE),
        footerStyle = list(fontWeight = "bold"),
        headerClass = "header"
      ),
      columns = list(
        age = colDef(
          minWidth = 200,
          name = "Age Group",
          footer = "Total (costs include Introduction cost)"
        ),
        cost_per_case = colDef(
          name = "Cost per Infection",
          format = colFormat(digits = 2),
          footer = function(values) {
            sprintf("$%s", format(round(total_cost_per_case, digits = 2), big.mark = ","))
          }
        ),
        total_incidence = colDef(
          name = "Total Incidence",
          footer = function(values) {
            result <- round(sum(values, na.rm = TRUE), digits = 0)
            sprintf("%s", format(result, big.mark = ",", nsmall = 0))
          }
        ),
        total_costs = colDef(
          name = "Total Costs",
          footer = function(values) {
            result <- round(sum(values, na.rm = TRUE) + input$cintro, digits = 0)
            sprintf("$%s", format(result, big.mark = ",", nsmall = 0))
          }
        )
      )
    )
  })

  # plot the model output: Incidence
  output$model_plot_incidence <- renderPlotly({
    plotInc(modelOutput() |> filter(time>=2022))
  })

  # plot the model output: Protected
  output$model_plot_protected <- renderPlotly({
    plotProt(modelOutput() |> filter(time>=2022))
  })

  # plot the model output: Treatment
  output$model_plot_treated <- renderPlotly({
    plotTr(modelOutput() |> filter(time>=2022))
  })
}

  # plot the model output: Treatment

shinyApp(ui, server)
