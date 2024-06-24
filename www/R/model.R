# model script

# Set the start and end time for the model simulation
times <- seq(0, 20, 1 / 365)

# set up a function to solve the equations
sirsv <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    betat <- beta * (1 + beta1 * cos(2 * pi * t))
    P <- S + I + R + V
    if (t > vacstart) {
      propvt <- propv/100
    } else {
      (propvt <- 0)
    }

    dS <- mu * P - betat * I / P * S * (1 - propvt) - mu * S + rho * R - v * propvt * S
    dI <- betat * I / P * S * (1 - propvt) - gamma * I - mu * I
    dR <- gamma * I - mu * R - rho * R
    dV <- v * propvt * S - mu * V
    dVacpop <- v * propvt * S
    output <- c(dS, dI, dR, dV, dVacpop)
    list(output)
  })
}