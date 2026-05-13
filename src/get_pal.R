get_pal <- function(x) {
  
  palettes <- list()
  
  palettes[["treatment"]] <- c("#FBA72A", "#5785C1")
  
  days <- c(
    "Before farrow",
    "Lactation",
    "D002",
    "D006",
    "D007",
    "D014",
    "D021",
    "D027",
    "D049",
    "D150"
  )
  
  palettes[["day"]] <- c(
    c("#899DA4", "#C93312"),
    RColorBrewer::brewer.pal(8, "Spectral")
  )
  
  names(palettes[["day"]]) <- days
  
  palettes[["hygiene"]] <- c(c("Optimal" = "#9986A5", "Basic" = "#79402E"))
  
  palettes[["feeders_per_pen"]] <- c(c("One" = "#90D4CC", "Two" = "#BD3027"))
  
  palettes[["pigs_per_pen"]] <- c(c("Seven" = "#85D4E3", "Fourteen" = "#F4B5BD"))
  
  palettes[["timepoint"]] <- c(
    "D0" = "#3B9AB2",
    "D14" = "#78B7C5",
    "D45" = "#EBCC2A",
    "Finisher" = "#E1AF00"
  )
  
  return(palettes[[x]])
  
}
