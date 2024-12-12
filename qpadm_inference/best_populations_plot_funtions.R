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

plot_best_2_pops <- function( population_rankings, all_ancestors, titlos = "" ){

  ## The following line extracts the 2 best population per simulation
  ## from a list of population rankings 
  best_2_populations <- melt(do.call(rbind, lapply(1:rep, function(x) {
    as.character(population_rankings[[x]][1:2, 1])
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
  barplot_of_best_2_pops <- barplot_of_best_2_pops + scale_x_discrete(drop = FALSE)
  
  return( barplot_of_best_2_pops )
}

# * Best Population Pair

## % of simulation that "population pair" was the best.

plot_best_pop_pair <- function( population_rankings, all_ancestors, titlos = "" ){

  ## The following line extracts the 2 best population per simulation
  ## from a list of population rankings 
  best_two_pops_2d <- as.data.frame(do.call(rbind, lapply(1:rep, function(x) {
    sort(as.character(population_rankings[[x]][1:2, 1]))
  })))
  names(best_two_pops_2d) <- c("best_1", "best_2")
  
  best_two_pops_2d$best_1 <- factor(best_two_pops_2d$best_1, levels = all_ancestors)
  best_two_pops_2d$best_2 <- factor(best_two_pops_2d$best_2, levels = all_ancestors)
  
  melted_best_two_pops_2d <- melt(table(best_two_pops_2d) / nrow(best_two_pops_2d)) # Divide by the number of simulations.
  triangle <- as.character(melted_best_two_pops_2d$best_1) < as.character(melted_best_two_pops_2d$best_2)
  
  heatmap_of_best_two_pops_2d <- ggplot(melted_best_two_pops_2d[triangle, ], aes(x = best_1, y = best_2, fill = value))
  heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + geom_tile(
    color = "white",
    lwd = 0.4,
    linetype = 1
  )
  heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + geom_text(
    aes(label = value),
    color = "white",
    size = 3
  )
  heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + labs(title = titlos)
  heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
  heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + scale_fill_gradient( low = "#56B4E9", high = "#E69F00")
  heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + scale_x_discrete(drop = FALSE)
  heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + scale_y_discrete(drop = FALSE)
  
  return( heatmap_of_best_two_pops_2d )
  
}


# * Accepted Models

## Number of models of "population" that are accepted divided by
## all the possible models that the population can take part in.

plot_accepted_models <- function( list_of_accepted_models, all_ancestors, titlos = "" ){

  average_accepted_models <- colMeans( do.call(rbind, list_of_accepted_models) )
  average_accepted_models_for_barplot <- data.frame(
    onomata = factor(names(average_accepted_models), levels = all_ancestors),
    values = as.numeric(sprintf("%.4f", average_accepted_models))
  )
  
  barplot_of_accepted_models <- ggplot(average_accepted_models_for_barplot, aes(x = onomata, y = as.numeric(values)))
  barplot_of_accepted_models <- barplot_of_accepted_models + geom_bar(stat = "identity", fill = "steelblue")
  barplot_of_accepted_models <- barplot_of_accepted_models + geom_text(
    aes(label = values),
    vjust = -0.3, size = 4
  )
  barplot_of_accepted_models <- barplot_of_accepted_models + labs(
    title = titlos
  )
  barplot_of_accepted_models <- barplot_of_accepted_models + theme(
    axis.title = element_blank()
  )
  barplot_of_accepted_models <- barplot_of_accepted_models + scale_x_discrete(drop = FALSE)

  return( barplot_of_accepted_models )
  
}

# * Accepted Models for Pair of Populations

## Number of models of "population pair" that are accepted divided by
## all the possible models that the population pair can take part in.

plot_accepted_models_2d <- funtion( list_of_accepted_models_2d, all_ancestors, titlos = "" ){

  melted_accepted_models_2d <- melt(Reduce("+", list_of_accepted_models_2d) / length(list_of_accepted_models_2d), na.rm = TRUE)
  melted_accepted_models_2d <- melted_accepted_models_2d[!is.nan(melted_accepted_models_2d$value), ]
  melted_accepted_models_2d$value  <- as.numeric(sprintf("%.4f", melted_accepted_models_2d$value))
  
  heatmap_of_accepted_models_plot_2d <- ggplot(melted_accepted_models_2d, aes( x = left_1, y = left_2, fill = value))
  heatmap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + geom_tile(
    color = "white",
    lwd = 0.4,
    linetype = 1
  )
  heatmap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + geom_text(
    aes(label = value),
    color = "white",
    size = 3
  )
  heatmap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + labs(title = "Scale 0.5")
  heatmapap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
  ## heatmap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + theme_minimal()
  heatmap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + scale_fill_gradient( low = "#56B4E9", high = "#E69F00")
  heatmap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + scale_y_discrete( drop = FALSE )
  heatmap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + scale_x_discrete( drop = FALSE )

  return( heatmap_of_accepted_models_2d )
  
}
