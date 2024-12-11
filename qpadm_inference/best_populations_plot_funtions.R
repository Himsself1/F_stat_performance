# * Packages

list_of_packages <- c(
  "ggplot2", "devtools",
  "argparse", "stringr",
  "Cairo", "tibble",
  "reshape", "dplyr"
)

for (i in list_of_packages) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i, dependencies = T)
  }
}

devtools::install_github("uqrmaie1/admixtools")
library(admixtools)

# * 2 Best Populations

## % of simulation that "population" was among the 2 best.

plot_best_2_pops <- function( population_comparisons, titlos = "" ){
  
  best_2_populations <- melt(do.call(rbind, lapply(1:rep, function(x) {
    as.character(population_comparisons[[x]][1:2, 1])
  })))
  
  best_2_for_barplot <- data.frame(
    onomata = factor(names(table(factor(best_2_populations[, 3], levels = all_ancestors))), levels = all_ancestors),
    values = as.vector(table(best_2_populations[, 3]))
  )
  
  barplot_of_best_2_pops <- ggplot(best_2_for_barplot, aes(x = onomata, y = as.numeric(values)))
  barplot_of_best_2_pops <- barplot_of_best_2_pops + geom_bar(stat = "identity", fill = "steelblue")
  barplot_of_best_2_pops <- barplot_of_best_2_pops + geom_text(
                                                       aes(label = values),
                                                       vjust = -0.3, size = 4
                                                     )
  barplot_of_best_2_pops <- barplot_of_best_2_pops + theme(
                                                       axis.title = element_blank()
                                                     )
barplot_of_best_2_pops <- barplot_of_best_2_pops + labs(
                                                     title = titlos
                                                   )
  }
