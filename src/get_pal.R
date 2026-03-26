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
  
  return(palettes[[x]])
  
}
