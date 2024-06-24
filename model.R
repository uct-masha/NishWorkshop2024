library(deSolve)
library(plotly)
library(tidyverse)
library(docstring) # Lets you use ?func for functions in this file

rates <- function(time, y, parms, contact) {
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
    Infectious <- c((In1 + It1), (In2 + It2), (In3 + It3))
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
    dVA2 <- cov1*a*R1 + cov1*a*S1 - ((1-ve1)*lambda[2] + a + mu)*VA2
    dR2 <- (1-cov1)*a*R1 + rr*In2 + rtr*Tr2 - a*R2 - mu*R2
    dCInc2 <- lambda[2]*S2 + (1-ve1)*lambda[2]*VA2
    dCTr2 <- rt*It2
    dCVax2 <- cov1*a*R1 + cov1*a*S1

    # Age group 3
    dS3 <- a*S2 - (lambda[3])*S3 - mu*S3
    dE3 <- a*E2 + lambda[3]*S3 + (1-ve1)*lambda[3]*VA3 - ((1-pt)*rs + pt*rs)*E3 - mu*E3
    dIn3 <- a*In2 + (1-pt)*rs*E3 - rr*In3 - mu*In3
    dIt3 <- a*It2 + pt*rs*E3 - rt*It3 - mu*It3
    dTr3 <- a*Tr2 + rt*It3 - rtr*Tr3 - mu*Tr3
    dVA3 <- - (1-ve1)*lambda[3]*VA3 +(1-cov2)*a*VA2 -mu*VA2
    dVB3 <- cov2*a*VA2 - mu*VB3
    dR3 <- a*R2+ rr*In3 + rtr*Tr3 - mu*R3
    dCInc3 <- lambda[3]*S3 + (1-ve1)*lambda[3]*VA3
    dCTr3 <- rt*It3
    dCVax3 <- cov2*a*VA2

    list(c(
      dS1, dE1, dIn1, dIt1, dTr1,             dR1, dCInc1, dCTr1,
      dS2, dE2, dIn2, dIt2, dTr2, dVA2,       dR2, dCInc2, dCTr2, dCVax2,
      dS3, dE3, dIn3, dIt3, dTr3, dVA3, dVB3, dR3, dCInc3, dCTr3, dCVax3
    ))
  })
}

runModel <- function(startyear=2010, endyear=2030, initialConditions, parameters, contact) {
  #' Run the model
  #' @param startyear The start year of the model
  #' @param endyear The end year of the model
  #' @param initialConditions The initial conditions of the model (flat vector)
  #' @param parameters The parameters of the model
  #' @param contact The contact matrix
  #' @return A tibble with the model results with `time`, `compartment`, and `population` columns
  timesteps <- seq(startyear, endyear, 1/365)
  mod <- ode(y = unlist(initialConditions),
             times = timesteps,
             func = rates,
             parms = parameters,
             contact = contact) |>
    as.data.frame() |>
    as_tibble() |>
    tidyr::pivot_longer(cols = !time, names_to = "compartment", values_to = "population") |>
    dplyr::mutate(compartment = factor(compartment, levels=names(initialConditions)))
  return(mod)
}

plotModel <- function(mod) {
  #' Plot the compartments of the model
  #'
  #' This function plots the compartments of the model
  #'
  #' @param mod The model to plot
  #' @return A plotly object
  plt <- ggplot2::ggplot(mod) +
    ggplot2::aes(x = time, y = population, color = compartment) +
    ggplot2::geom_line() +
    labs(title = "Compartmental Model",
         color = "Compartment",
         x = "Time",
         y = "Population")
  ggplotly(plt)
}

plotPop <- function(mod) {
  plt <- mod |>
    filter(compartment |> str_starts("C", negate = T)) |>
    mutate(age=str_extract(compartment,pattern="\\d+$") |> as.numeric()) |>
    summarise(pop=sum(population), .by=c(time,age)) |>
    ggplot(aes(x=time, y=pop, color=factor(age))) +
    geom_line() +
    labs(title="Population by age group",
         x="Time",
         y="Population",
         color="Age group")
  ggplotly(plt)
}

plotInc <- function(mod, byAge=TRUE) {
  tbIncAges <- mod |>
    mutate(age=str_extract(compartment,pattern="\\d+$") |> as.numeric(),
           Year=year(date_decimal(time)),
           Month=month(date_decimal(time))) |>
    filter(compartment |> str_starts("CInc")) |>
    summarise(Incidence=last(population) - first(population),
              .by=c(Year, age))
  tbIncAges |>
    ggplot(aes(x=Year, y=Incidence, fill=factor(age))) +
    geom_col(position="dodge") +
    labs(title="Incidence by age group",
         x="Year",
         y="Incidence",
         fill="Age group")
}

plotTr <- function(mod, byAge=TRUE) {
  tbIncAges <- mod |>
    mutate(age=str_extract(compartment,pattern="\\d+$") |> as.numeric(),
           Year=year(date_decimal(time)),
           Month=month(date_decimal(time))) |>
    filter(compartment |> str_starts("CTr")) |>
    summarise(Incidence=last(population) - first(population),
              .by=c(Year, age))
  tbIncAges |>
    ggplot(aes(x=Year, y=Incidence, fill=factor(age))) +
    geom_col(position="dodge") +
    labs(title="Treatment by age group",
         x="Year",
         y="Treatment",
         fill="Age group")
}

makeInitialConditions <- function(
    S1=100,  E1=1, In1=1, It1=1, Tr1=1,               R1=0,
    S2=110,  E2=1, In2=1, It2=1, Tr2=1, VA2=0,        R2=0,
    S3=100, E3=1, In3=1, It3=1, Tr3=1, VA3=0, VB3=0,  R3=900) {
  list(S1=S1, E1=E1, In1=In1, It1=It1, Tr1=Tr1,                  R1=R1, CInc1=0, CTr1=0,
       S2=S2, E2=E2, In2=In2, It2=It2, Tr2=Tr2, VA2=VA2,         R2=R2, CInc2=0, CTr2=0, CVax2=0,
       S3=S3, E3=E3, In3=In3, It3=It3, Tr3=Tr3, VA3=VA3, VB3=VB3,R3=R3, CInc3=0, CTr3=0, CVax3=0)
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
  c(sum(x[comp1]),
    sum(x[comp2]),
    sum(x[comp3]))
}

# survey_zim <- socialmixr::get_survey("https://doi.org/10.5281/zenodo.3886638")
# cm <- socialmixr::contact_matrix(survey_zim, countries="Zimbabwe",
#                                  age.limits = c(0,1,2),
#                                  sample.participants = T,
#                                  symmetric = T,
#                                  return.demography = T)
#
# Mp=matrix(0,nrow=3,ncol=3)
# for (i in 1:3) {
#   for (j in 1:3) {
#     Ni <- cm$demography$population[i]
#     Nj <- cm$demography$population[j]
#     Npi <- pop[i]
#     Npj <- pop[j]
#     Np <- sum(pop)
#     Mij <- cm$matrix[i,j]
#     sigmaDenom <- sum(pop*cm$matrix[i,]) # TODO: Fix this denom from Arregui M3
#     Mp[i,j] <- Mij * (Npj/Nj) * Np/sigmaDenom
#   }
# }
# cm

makeContactMatrix <- function(initialConditions, total_contacts=matrix(c(11,  5,   100,
                                                                         5,   18,  400,
                                                                         100, 400, 2e4), nrow=3)*365) {
  pop <- sum(as.double(initialConditions))
  ageProps <- with(initialConditions,
                   c(S1+E1+In1+It1+Tr1+        R1,
                     S2+E2+In2+It2+Tr2+VA2+    R2,
                     S3+E3+In3+It3+Tr3+VA3+VB3+R3)/pop
  )
  contact <- total_contacts/pop
  return(contact)
}

makeParameters <- function(mu=0.01,
                           beta=0.1,
                           cov1=0.65,
                           cov2=0.5,
                           pt=1,
                           rs=1,
                           rr=1,
                           ve1=1,
                           rt=1,
                           rtr=1,
                           a=1,  # a=ageing=1yr
                           timevax=2025) {
  c(mu=mu,
    beta=beta,
    cov1=cov1,
    cov2=cov2,
    pt=pt,
    rs=rs,
    rr=rr,
    ve1=ve1,
    rt=rt,
    rtr=rtr,
    a=a,
    timevax=timevax
  )
}

initialConditions <- makeInitialConditions()
total_contacts <- matrix(c(20, 2, 1,
                           2, 30, 4,
                           1, 4, 80), nrow=3)*365

contact <- makeContactMatrix(initialConditions = initialConditions, total_contacts = total_contacts)

mod <- runModel(initialConditions=initialConditions, parameters=makeParameters(), contact = contact)

# plotModel(mod)
plotInc(mod |> filter(time>=2024))
