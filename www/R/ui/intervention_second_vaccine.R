div(
  class = "interventions",
  fluidRow(
    column(
      width = 8,
      materialSwitch(
        inputId = "prim_switch2", label = "Vaccine 2", value = FALSE,
        status = "info", right = TRUE, inline = FALSE, width = NULL
      )
    ),
    column(
      width = 4,
      conditionalPanel(
        condition = "input.prim_switch2",
        dropdownButton(
          label = "Settings", circle = FALSE, status = "primary", size = "sm", icon = icon("gear"),
          width = "800px", margin = "15px", right = TRUE,
          layout_column_wrap(
            width = 1/2,
            sliderInput(inputId = "yr_start2", label = "Year to Start", value = 2026, min = 2025, max = 2030, sep = ""),
            sliderInput(inputId = "cov2", label = "Vaccine coverage", value = 60, min = 1, max = 100),
            sliderInput(inputId = "ytrc2", label = "Years to reach coverage:", value = 2, min = 1, max = 10),
            selectInput(inputId = "delivery2", label = "Method of delivery", choices = c("Health facility", "Outreach site"))
          )
        )
      )
    )
  ),
  tags$small("This is a detailed description for the second vaccine.")
)