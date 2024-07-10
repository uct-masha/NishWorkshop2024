source("model.R")
source("packages.R")

initialConditions <- rescalePop(unlist(makeInitialConditions()), 50000000)
contact <- makeContactMatrix(initialConditions = initialConditions)

mod <- runModel(startyear=2000, endyear=2040,
                initialConditions=initialConditions,
                parameters=makeParameters(infectiousWhileTreated = 0),
                contact = contact)

tableCosts(mod = mod, cintro = 500000, cvacc = 2, cdel = 3, ctrt = 1, discountRate = 0.03, discountYear = 2025)

# plotModel(mod, Variables = c("Inc","Tr","Vax"))
# plotModel(mod, Variables = c("Inc"), Ages = 2:4)
# plotModel(mod, Variables = c("Pop"))
# plotModel(mod, Compartments = c("S","E","In","It","Tr","VA","VB","R","CInc","CTr","CVax"))
# plotModel(mod, Compartments = c("VA", "VB"), Ages = 2:4)
