library(grow96)

#glucose <- c(0.2, 0.15, 0.1, 0.05, 0.01, 0)
glucose <- c(0.2)
M9_MGSO4 <- c(100)
NH4CL <- c(0, 20, 50, 100, 150, 200)
volume <- c(20, 60, 100, 140, 180, 100)

overnight_technical <- outer(paste0("O", 1:3, "_"), paste0("T", 1:3), FUN = "paste0")
attributes(overnight_technical) <- NULL
overnight_technical <- c(overnight_technical, NA)

makeSpec_fullFact(plateName = "AT_oxy_gradient",
                  rowName = "volume",
                  columnName = "overnight_technical",
                  rows = volume,
                  columns = overnight_technical,
                  border = "EMPTY",
                  blankColumns = 11,
                  replicates = 1,
                  constants = c("glucose" = 0.2, "M9_MGSO4" = 100),
                  specPath = "./specs")
