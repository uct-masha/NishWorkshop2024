
rates <- function(time, y, parms) {
  #' Model rates function
  #'
  #' This function calculates the rates of change for the compartments in the model.
  #' @param time The current time
  #' @param y The current state of the model
  #' @param parms The parameters of the model
  #' @param contact The contact matrix
  #' @return A list of the rates of change for each compartment
  popVec <- getPopVec(y)
  pop <- sum(popVec)
  with(as.list(c(y, parms)), {
    if (parms['infectiousWhileTreated']!=0) {
      Infectious <- c((In1 + It1 + Tr1), (In2 + It2 + Tr2), (In3 + It3 + Tr3), (In4 + It4 + Tr4))
    } else {
      Infectious <- c((In1 + It1), (In2 + It2), (In3 + It3), (In4 + It4))
    }
    lambda <- as.double(beta*contact%*%Infectious/popVec)
    if (time <= timevax) {cov1 <- cov2 <- 0}

    # Age group 1
    dS1 <- mu*pop - (lambda[1] + a+mu)*S1
    dE1 <- lambda[1]*S1 - ((1-pt)*rs + pt*rs +a+mu)*E1
    dIn1 <- (1-pt)*rs*E1 - (rr+a+mu)*In1
    dIt1 <- pt*rs*E1 - (rt+a+mu)*It1
    dTr1 <- rt*It1 - (rtr+a+mu)*Tr1
    dR1 <- rr*In1 + rtr*Tr1 - (a+mu)*R1
    dCInc1 <- lambda[1]*S1
    dCTr1 <- rt*It1

    # Age group 2
    dS2 <- (1-cov1)*a*S1 -(lambda[2] + a +mu)*S2
    dE2 <- a*E1 + lambda[2]*S2 + (1-ve1)*lambda[2]*VA2 - ((1-pt)*rs + pt*rs +a +mu)*E2
    dIn2 <- a*In1 + (1-pt)*rs*E2 - (rr+a +mu)*In2
    dIt2 <- a*It1 + pt*rs*E2 - (rt+a +mu)*It2
    dTr2 <- a*Tr1 + rt*It2 - (rtr+a + mu)*Tr2
    dVA2 <- cov1*a*S1 - ((1-ve1)*lambda[2] + a + mu)*VA2
    dR2 <- a*R1 + rr*In2 + rtr*Tr2 - a*R2 - mu*R2
    dCInc2 <- lambda[2]*S2 + (1-ve1)*lambda[2]*VA2
    dCTr2 <- rt*It2
    dCVax2 <- cov1*a*S1

    # Age group 3
    dS3 <- (1-cov2)*a*S2 - (lambda[3])*S3 - (mu+a2)*S3 + tau1*VA3
    dE3 <- a*E2 + lambda[3]*S3 + (1-ve1)*lambda[3]*VA3 - ((1-pt)*rs + pt*rs)*E3 - (mu+a2)*E3
    dIn3 <- a*In2 + (1-pt)*rs*E3 - rr*In3 - (mu+a2)*In3
    dIt3 <- a*It2 + pt*rs*E3 - rt*It3 - (mu+a2)*It3
    dTr3 <- a*Tr2 + rt*It3 - rtr*Tr3 - (mu+a2)*Tr3
    dVA3 <- cov2*a*S2 + (1-cov2)*a*VA2 - (1-ve1)*lambda[3]*VA3 -(mu + a2 + tau1)*VA3
    dVB3 <- cov2*a*VA2 - (mu+a2)*VB3
    dR3 <- a*R2+ rr*In3 + rtr*Tr3 - (mu+a2)*R3
    dCInc3 <- lambda[3]*S3 + (1-ve1)*lambda[3]*VA3
    dCTr3 <- rt*It3
    dCVax3 <- cov2*a*S2 +cov2*a*VA2

    # Age group 4
    dS4 <- a2*S3 - (lambda[4])*S4 - mu*S4 + tau1*VA4
    dE4 <- a2*E3 + lambda[4]*S4 + (1-ve1)*lambda[4]*VA4 - ((1-pt)*rs + pt*rs)*E4 - mu*E4
    dIn4 <- a2*In3 + (1-pt)*rs*E4 - rr*In4 - mu*In4
    dIt4 <- a2*It3 + pt*rs*E4 - rt*It4 - mu*It4
    dTr4 <- a2*Tr3 + rt*It4 - rtr*Tr4 - mu*Tr4
    dVA4 <- a2*VA3 - (1-ve1)*lambda[4]*VA4 -(mu + tau1)*VA4
    dVB4 <- a2*VB3 - mu*VB4
    dR4 <- a2*R3 + rr*In4 + rtr*Tr4 - mu*R4
    dCInc4 <- lambda[4]*S4 + (1-ve1)*lambda[4]*VA4
    dCTr4 <- rt*It4

    list(c(
      dS1, dE1, dIn1, dIt1, dTr1,             dR1, dCInc1, dCTr1,
      dS2, dE2, dIn2, dIt2, dTr2, dVA2,       dR2, dCInc2, dCTr2, dCVax2,
      dS3, dE3, dIn3, dIt3, dTr3, dVA3, dVB3, dR3, dCInc3, dCTr3, dCVax3,
      dS4, dE4, dIn4, dIt4, dTr4, dVA4, dVB4, dR4, dCInc4, dCTr4
    ))
  })
}

postProc <- function(mod) {
  #' Post-process the model
  #' @param mod The model to post-process
  #' @return A tibble with the post-processed results
  popage <- mod |> #create average population by age and year for denominator in plots
    filter(compartment |> str_starts("C", negate = T)) |>
    mutate(age=str_extract(compartment,pattern="\\d+$") |> as.numeric()) |>
    summarise(pop=sum(population), .by=c(time,age)) |>
    mutate(Year = floor(time)) |>
    summarise(popyr=((last(pop) + first(pop))/2),
              .by=c(Year, age))

}

ageGroupNames <- function(){"0-1, 1-2, 2-5, 5+" |> str_split(", ") |> unlist()}

runModel <- function(startyear=2000, endyear=2040, initialConditions, parameters, contact) {
  #' Run the model
  #' @param startyear The start year of the model
  #' @param endyear The end year of the model
  #' @param initialConditions The initial conditions of the model (flat vector)
  #' @param parameters The parameters of the model
  #' @param contact The contact matrix
  #' @return a list with the daily and annual result tibbles. The daily tibble has the following columns: time, compartment, population, totpop. The annual tibble has the following columns: Year, compartment, population
  timesteps <- seq(startyear, endyear, 1/365)
  # The parameters need to be a list and they need to include contact
  parms = as.list(parameters)
  parms$contact = contact

  moRaw <- dde(y = unlist(initialConditions),
             times = timesteps,
             func = rates,
             parms = parms)

  # Put into long format tibble with Pop "compartment" used for per capita calculations
  moDaily <- moRaw |>
    as.data.frame() |>
    as_tibble() |>
    pivot_longer(cols = !time, names_to = "compartment", values_to = "population") |>
    # Split out compartment into compartment and age
    separate_wider_regex(compartment, c(compartment="[a-zA-Z]+", AgeGroup="\\d+")) |>
    # Include Year and improve format of AgeGroup and compartment columns
    mutate(Year=floor(time),
           AgeGroup = as_factor(ageGroupNames()[as.numeric(AgeGroup)]),
           compartment = as_factor(compartment)) %>%
    bind_rows(., # Include rows for "Pop" variable containing total alive population per age group
              . |> summarise(compartment="Pop", population=sum(population[!str_starts(compartment,"C")]), .by=c(time, AgeGroup, Year)))
    # Include totpop column with annual age-specific population
    # reframe(time=time,
    #         compartment=compartment,
    #         population=population,
    #         totpop=sum(population[compartment |> str_starts("C", negate = T)])/365,
    #         .by=c(Year,AgeGroup))

  # Calculate annual values: incidence, treatment, vaccination
  moAnnual <- moDaily |>
    filter(str_starts(compartment, "C")|compartment=="Pop") |>
    reframe(population = if_else(compartment=="Pop",
                                   mean(population),
                                   diff(range(population))),
              .by = c(Year, compartment, AgeGroup)) |>
    distinct() |>
    mutate(compartment = str_remove(compartment, "^C")) |>
    pivot_wider(names_from = compartment, values_from = population)

  return(list(
    moDaily = moDaily,
    moAnnual = moAnnual
  ))
}

plotModel <- function(mod,
                      Variables=NULL,
                      Compartments=NULL,
                      Ages=1:4,
                      TimeRange=c(2000,2050),
                      Per1KPop=F,
                      UsePlotly=F) {
  #' Plot model outputs
  #'
  #' This function plots the compartments of the model generically.
  #' Exactly one of `Variables` or `Compartments` must be given.
  #'
  #' @param mod The model to plot (list with daily and annual tibbles)
  #' @param Variables The variables to plot. One of `c("Inc", "Tr", "Vax", "Pop")`
  #' @param Compartments The compartments to plot. One of `"S", "It", ..., "Pop"`.
  #' @param Ages The ages to show in the plot. Must be a subset of 1:4 of length>=1.
  #' @param TimeRange The years to plot between
  #' @param Per1KPop Whether to plot per 1000 population (if `Variables` is given)
  #' @param UsePlotly Whether to use plotly for the plot
  #' @return A plotly object

  # Validate inputs
  if (missing(Variables) & missing(Compartments)) {
    stop("One of Variables or Compartments must be given")
  }
  if (!missing(Variables) & !missing(Compartments)) {
    stop("Only one of Variables or Compartments can be given")
  }
  if ("Pop" %in% Variables & length(Variables) > 1) {
    stop("Pop cannot be plotted with other variables. Plot it alone or only use other variables")
  }
  if ("Pop" %in% Compartments & length(Compartments) > 1) {
    stop("Pop cannot be plotted with other compartments. Plot it alone or only use other compartments")
  }
  if (!all(Ages %in% 1:4)) {
    stop("Ages must be a subset of 1:4")
  }
  if (!missing(Variables)) {
    if (!all(Variables %in% c("Inc", "Tr", "Vax", "Pop"))) {
      stop("Variables must be a subset of c('Inc', 'Tr', 'Vax', 'Pop')")
    }
  } else {
    if (!all(Compartments %in% c("S", "E", "In", "It", "Tr", "VA", "VB", "R", "CInc", "CTr", "CVax"))) {
      stop("Compartments must be a subset of c('S', 'E', 'In', 'It', 'Tr', 'VA', 'VB', 'R', 'CInc', 'CTr', 'CVax')")
    }
  }

  shouldMakeAnnualPlot <- !missing(Variables)
  # Make the plots as required
  if (shouldMakeAnnualPlot) {
    if ("Pop" %in% Variables) {
      # User is asking to plot population only
      plt <- mod$moAnnual |>
        filter(as.integer(AgeGroup) %in% Ages,
               Year >= first(TimeRange) & Year <= last(TimeRange)) |>
        ggplot(aes(x = Year, y = Pop, color = factor(AgeGroup))) +
        geom_line() +
        labs(
          title = paste0(popTxt, " by age group"),
          x = "Year",
          y = popTxt,
          color = "Age group")
    } else {
      # User is asking to plot non-population variables (Inc, Tr, Vax)
      popTxt <- if(Per1KPop) {"Rate per 1000 population"}else{"Population"}
      plt <- mod$moAnnual |>
        filter(as.integer(AgeGroup) %in% Ages,
               Year >= first(TimeRange) & Year <= last(TimeRange)) |>
        mutate(popScale = ifelse(Per1KPop,1000/Pop,1)) |>
        pivot_longer(cols = all_of(Variables), names_to = "variable", values_to = "pop") |>
        mutate(
          variableName = case_when(
            variable == "Inc" ~ "Incidence",
            variable == "Tr" ~ "Treated",
            variable == "Vax" ~ "Vaccinated"
          )
        ) |>
        select(Year, AgeGroup, variableName, pop) |>
        ggplot(aes(x = Year,
                   y = pop,
                   group = AgeGroup,
                   fill = AgeGroup)) +
        geom_col(position = "dodge") +
        facet_wrap(~variableName, scales = "free_y") +
        scale_y_continuous(labels = scales::comma) +
        labs(
          title = paste0(popTxt, " by age group"),
          x = "Year",
          y = popTxt,
          color = "Age group")
    }
  } else {
    # Make the daily plot
    plt <- mod$moDaily |>
      filter(as.integer(AgeGroup) %in% Ages,
             compartment %in% Compartments,
             Year >= first(TimeRange) & Year <= last(TimeRange)) |>
      mutate(compartment=as_factor(compartment)) |>
      ggplot(aes(x = time,
                 y = population,
                 group = AgeGroup,
                 color = AgeGroup)) +
      geom_line() +
      facet_wrap(~compartment, scales = "free_y") +
      scale_y_continuous(labels = scales::comma) +
      labs(
        title = "Population by age group",
        x = "Time",
        y = "Population",
        color = "Age group"
      )
  }
  if (UsePlotly) {
    plotly::ggplotly(plt)
  } else {
    plt
  }
}

# plotModel(mod, Variables = c("Inc","Tr","Vax"))
# plotModel(mod, Variables = c("Inc"), Ages = 2:4)
# plotModel(mod, Variables = c("Pop"))
# plotModel(mod, Compartments = c("S","E","In","It","Tr","VA","VB","R","CInc","CTr","CVax"))
# plotModel(mod, Compartments = c("VA", "VB"), Ages = 2:4)

getCosts <- function(mod, cintro, cvacc, cdel, ctrt) {
  #' Get the costs of the model
  #' @param mod The model to get costs from
  #' @return A tibble with the costs

  ageGroupsVaccinated <- mod$moAnnual |>
    filter(Vax>0) |>
    pull(AgeGroup) |>
    unique() |>
    as.integer()

  mod$moAnnual |>
    mutate(
      Doses = if_else(is.na(Vax), 0, Vax),
      CostIntro = if_else(as.integer(AgeGroup) %in% ageGroupsVaccinated, cintro, 0),
      CostVax = cvacc * Doses,
      CostDel = cdel * Doses,
      CostTr = ctrt * Tr,
      CostTot = Doses + CostDel + CostTr
    ) |>
    select(Year, AgeGroup, Doses:CostTot) |>
    pivot_longer(cols=!c(Year,AgeGroup), names_to="variable", values_to="Value")
}

getCostsEpi <- function(mod, cintro, cvacc, cdel, ctrt, discountRate=0.03, discountYear=2025) {
  getCosts(mod, cintro, cvacc, cdel, ctrt) |>
    left_join(mod$moAnnual, by = join_by(Year, AgeGroup)) |>
    pivot_wider(names_from = "variable", values_from = "Value") |>
    mutate(discounting_factor = (1 + discountRate) ^ (discountYear - Year),
           TotalCostsDiscounted = CostTot * discounting_factor) |>
    summarise(TotalCosts = sum(TotalCostsDiscounted),
              TotalIncidence = sum(Inc),
              CostPerCase = TotalCosts / TotalIncidence,
              .by = AgeGroup)
}

tableCosts <- function(mod, cintro, cvacc, cdel, ctrt, discountRate=0.03, discountYear=2025) {
  getCostsEpi(mod, cintro, cvacc, cdel, ctrt, discountRate, discountYear) |>
    reactable(defaultColDef = colDef(
                format = colFormat(digits = 0, separators = TRUE),
                footerStyle = list(fontWeight = "bold")
              ),
              columns = list(
                AgeGroup = colDef(
                  minWidth = 200,
                  name = "Age Group",
                  footer = "Total (costs include Introduction cost)"
                ),
                CostPerCase = colDef(
                  name = "Cost per Incidence",
                  format = colFormat(digits = 2),
                  footer = function(values) {
                    result <- round(sum(values, na.rm = TRUE), digits = 2)
                    sprintf("$%s", format(result, big.mark = ","))
                  }
                ),
                TotalIncidence = colDef(
                  name = "Total Incidence",
                  footer = function(values) {
                    result <- round(sum(values, na.rm = TRUE), digits = 0)
                    sprintf("$%s", format(result, big.mark = ",", nsmall = 0))
                  }
                ),
                TotalCosts = colDef(
                  name = "Total Costs",
                  footer = function(values) {
                    result <- round(sum(values, na.rm = TRUE) + cintro, digits = 0)
                    sprintf("$%s", format(result, big.mark = ",", nsmall = 0))
                  }
                )
              )
    )

}

# tableCosts(mod = mod, cintro = 500000, cvacc = 2, cdel = 3, ctrt = 1, discountRate = 0.03, discountYear = 2025)

# getCosts(mod, cintro=500000, cvacc=2, cdel=3, ctrt=1)

costing <- function(mod, byAge=TRUE) {
  tbOutAges <- mod |>
    mutate(age=str_extract(compartment,pattern="\\d+$") |> as.numeric(),
           Year=year(date_decimal(time)),
           Month=month(date_decimal(time))) |>
    filter(compartment |> str_starts("C"), Year >= makeParameters()["timevax"]) |>
    summarise(Value=last(population) - first(population),
              .by=c(Year, age, compartment)) |>
   mutate(variable = case_when(
     str_starts(compartment, "CInc") ~ "Inc",
     str_starts(compartment,"CTr") ~ "Tr",
     str_starts(compartment, "CVax") ~ "Vax"
    )) |>
    select(-compartment) |>
    pivot_wider(names_from=variable, values_from=Value) |>
    mutate(
      CostIntro = makeParameters()["cintro"],
      CostVax = Vax * makeParameters()["cvacc"],
      CostDel = Vax * makeParameters()["cdel"],
      CostTr = Tr * makeParameters()["ctrt"],
      CostTot = if_else(is.na(CostVax), 0, CostVax) + if_else(is.na(CostDel), 0, CostDel) + if_else(is.na(CostTr), 0, CostTr)
    ) |>
    pivot_longer(cols=c(Inc, Tr, Vax, CostIntro, CostDel, CostVax, CostTr, CostTot), names_to="variable", values_to="Value")
}

makeInitialConditions <- function(
    S1=96,  E1=1, In1=5, It1=1, Tr1=1,               R1=0,
    S2=98,  E2=1, In2=3, It2=1, Tr2=1, VA2=0,        R2=0,
    S3=105, E3=1, In3=1, It3=1, Tr3=1, VA3=0, VB3=0, R3=200,
    S4=500, E4=1, In4=1, It4=1, Tr4=1, VA4=0, VB4=0, R4=5000) {
  list(S1=S1, E1=E1, In1=In1, It1=It1, Tr1=Tr1,                  R1=R1, CInc1=0, CTr1=0,
       S2=S2, E2=E2, In2=In2, It2=It2, Tr2=Tr2, VA2=VA2,         R2=R2, CInc2=0, CTr2=0, CVax2=0,
       S3=S3, E3=E3, In3=In3, It3=It3, Tr3=Tr3, VA3=VA3, VB3=VB3,R3=R3, CInc3=0, CTr3=0, CVax3=0,
       S4=S4, E4=E4, In4=In4, It4=It4, Tr4=Tr4, VA4=VA4, VB4=VB4,R4=R4, CInc4=0, CTr4=0)
}

rescalePop <- function(popListCurrent, popTotal) {
  popTotalCurrent <- sum(popListCurrent)
  lapply(popListCurrent, function(x) (x*popTotal)/popTotalCurrent)
}



# mod |>
#   filter(time==max(time)) |> pull(population, name = compartment) |>
#   saveRDS("ic.rds")
# initialConditions <- as.list(readRDS("ic.RDS"))

getPopVec <- function(x) {
  nx <- names(x)
  notC <- str_starts(nx, "C", negate = T)
  comp1 <- notC & str_ends(nx, "1")
  comp2 <- notC & str_ends(nx, "2")
  comp3 <- notC & str_ends(nx, "3")
  comp4 <- notC & str_ends(nx, "4")

  c(sum(x[comp1]),
    sum(x[comp2]),
    sum(x[comp3]),
    sum(x[comp4]))
}

makeContactMatrix <- function(initialConditions) {
  #' Roughly put in an estimate of contacts per day
  #' Keep in mind the size of the age groups
  #' In our model: 0-1, 1-2, 2-5, 5+
  #' @param initialConditions The initial conditions of the model
  #' @return A contact matrix which preserves reciprocity
  pops <- purrr::map_dbl(1:4, ~sum(as.numeric(initialConditions[stringr::str_detect(names(initialConditions), as.character(.x))])))
  contact <- matrix(c(
    0.48, 0.30, 0.36, 1.6,  # 0-1 year old has how many contacts with other age group?
    0.30, 0.84, 0.41, 1.2,
    0.36, 0.41, 1.61, 2.2,
    0.01, 0.02, 0.3, 1.9
  ), byrow = T, ncol=4)
  # Naive way to ensure reciprocity - one should investigate the terms in the sum
  contact <- (contact*pops + t(contact*pops))/2/pops
  contact
}

makeParameters <- function(mub= 1/45,
                           mu=1/65,
                           beta=61,
                           cov1=0.0,
                           cov2=0.0,
                           pt=0.7,
                           rs=365/10,
                           rr=365/30,
                           ve1=0.85, #Vaccine efficacy dose A; dose B = 100%
                           tau1 = 0, #No waning as vaccine offers lifelong protection
                           rt=365/2,
                           rtr=365/7,
                           a=1,  # a=ageing=1yr
                           a2 = 1/3, #a2 = ageing over 3 years
                           timevax=2025,
                           discount = 0.03,
                           cvacc = 1.5,
                           cdel = 1,
                           ctrt = 0.5,
                           cintro = 500000,
                           infectiousWhileTreated = 0) {
  c(mub=mub,
    mu=mu,
    beta=beta,
    cov1=cov1,
    cov2=cov2,
    pt=pt,
    rs=rs,
    rr=rr,
    ve1=ve1,
    tau1=tau1,
    rt=rt,
    rtr=rtr,
    a=a,
    a2=a2,
    timevax=timevax,
    discount = 0.03,
    cvacc = cvacc,
    cdel = cdel,
    ctrt = ctrt,
    cintro = cintro,
    infectiousWhileTreated = infectiousWhileTreated
  )
}

