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

plot_best_2_pops <- function( population_rankings, good_models, all_ancestors, titlos = "" ){

  filtered_rankings <- population_rankings[good_models > 0]
  ## The following line extracts the 2 best population per simulation
  ## from a list of population rankings 
  best_2_populations <- melt(do.call(rbind, lapply(1:length(filtered_rankings), function(x) {
    as.character(filtered_rankings[[x]][1:2, 1])
  })))
  
  best_2_for_barplot <- data.frame(
    onomata = factor(names(table(factor(best_2_populations[, 3], levels = all_ancestors))), levels = all_ancestors),
    values = as.vector(table(factor(best_2_populations[, 3], levels = all_ancestors))/length(good_models))
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
  barplot_of_best_2_pops <- barplot_of_best_2_pops + scale_y_continuous(limits = c(0.0, 1.0))
  return( barplot_of_best_2_pops )
}

# * Best Population Pair

## % of simulation that "population pair" was the best.

plot_best_pop_pair <- function( population_rankings, good_models, all_ancestors, titlos = "" ){

  filtered_rankings <- population_rankings[good_models > 0]
  ## The following line extracts the 2 best population per simulation
  ## from a list of population rankings 
  best_two_pops_2d <- as.data.frame(do.call(rbind, lapply(1:length(filtered_rankings), function(x) {
    sort(as.character(filtered_rankings[[x]][1:2, 1]))
  })))
  names(best_two_pops_2d) <- c("best_1", "best_2")
  
  best_two_pops_2d$best_1 <- factor(best_two_pops_2d$best_1, levels = all_ancestors)
  best_two_pops_2d$best_2 <- factor(best_two_pops_2d$best_2, levels = all_ancestors)
  
  melted_best_two_pops_2d <- melt(table(best_two_pops_2d) / length(good_models)) # Divide by the number of simulations.
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
  heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + scale_fill_gradient2( low = "#56B4E9", high = "#E69F00", midpoint = 0.5, limits = c(0.0, 1.0))
  heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + scale_x_discrete(drop = FALSE)
  heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + scale_y_discrete(drop = FALSE)
  
  return( heatmap_of_best_two_pops_2d )
  
}

# * Accepted Models

## Number of models of "population" that are accepted divided by
## all the possible models that the population can take part in.

plot_accepted_models <- function( list_of_accepted_models, good_models, all_ancestors, titlos = "" ){

  filtered_accepted_models <- list_of_accepted_models[good_models > 0]
  average_accepted_models <- colSums( do.call(rbind, filtered_accepted_models) )/length(good_models)
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
  barplot_of_accepted_models <- barplot_of_accepted_models + scale_y_continuous(limits = c(0.0, 1.0))
  return( barplot_of_accepted_models )
  
}

# * Accepted Models for Pair of Populations

## Number of models of "population pair" that are accepted divided by
## all the possible models that the population pair can take part in.

plot_accepted_models_2d <- function( list_of_accepted_models_2d, good_models, all_ancestors, titlos = "" ){

  filtered_accepted_models_2d <- list_of_accepted_models_2d[good_models > 0]
  ## melted_accepted_models_2d <- melt(Reduce("+", filtered_accepted_models_2d) / length(filtered_accepted_models_2d), na.rm = TRUE)
  melted_accepted_models_2d <- melt(Reduce("+", filtered_accepted_models_2d) / length(good_models), na.rm = TRUE)
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
  heatmap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + labs(title = titlos)
  ## heatmapap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + theme(
  ##   axis.title.x = element_blank(),
  ##   axis.title.y = element_blank()
  ## )
  heatmap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + theme(
      axis.title = element_blank()
    )
  ## heatmap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + theme_minimal()
  heatmap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + scale_fill_gradient2( low = "#56B4E9", high = "#E69F00", midpoint = 0.5, limits = c(0.0, 1.0) )
  heatmap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + scale_y_discrete( drop = FALSE )
  heatmap_of_accepted_models_plot_2d <- heatmap_of_accepted_models_plot_2d + scale_x_discrete( drop = FALSE )

  return( heatmap_of_accepted_models_plot_2d )
  
}

# * Specificity Plots

plot_specificity <- function( list_of_specificity, good_models, all_ancestors, titlos = "" ){
  
  filtered_specificity <- list_of_specificity[good_models > 0]
  average_specificity <- colSums(do.call(rbind, filtered_specificity))/length(good_models)
  average_specificity_for_barplot <- data.frame(
    onomata = factor(names(average_specificity), levels = all_ancestors),
    values = as.numeric(sprintf("%.4f", average_specificity))
  )
  
  barplot_of_specificity <- ggplot(average_specificity_for_barplot, aes(x = onomata, y = as.numeric(values)))
  barplot_of_specificity <- barplot_of_specificity + geom_bar(stat = "identity", fill = "steelblue")
  barplot_of_specificity <- barplot_of_specificity + geom_text(
    aes(label = values),
    vjust = -0.3, size = 4
  )
  barplot_of_specificity <- barplot_of_specificity + labs(
    title = titlos
  )
  barplot_of_specificity <- barplot_of_specificity + theme(
    axis.title = element_blank()
  )
  barplot_of_specificity <- barplot_of_specificity + scale_x_discrete(drop = FALSE)
  barplot_of_specificity <- barplot_of_specificity + scale_y_continuous( limits(0.0, 1.0) )
  return( barplot_of_specificity )
  
}

# * Specificity Pair Plots

plot_specificity_2d <- function( list_of_specificity_2d, good_models, all_ancestors, titlos = "" ){

  filtered_specificity_2d <- list_of_specificity_2d[good_models > 0]
  temp_reduce <- Reduce("+", filtered_specificity_2d) / length(filtered_specificity_2d)
  ## melted_specificity_2d <- melt(temp_reduce / (length(filtered_specificity_2d) * upper.tri(temp_reduce)), na.rm = TRUE )
  melted_specificity_2d <- melt(temp_reduce / (length(good_models) * upper.tri(temp_reduce)), na.rm = TRUE )
  melted_specificity_2d <- melted_specificity_2d[!is.nan(melted_specificity_2d$value), ]
  melted_specificity_2d$value  <- as.numeric(sprintf("%.4f", melted_specificity_2d$value))
  
  heatmap_of_specificity_plot_2d <- ggplot(melted_specificity_2d, aes( x = left_1, y = left_2, fill = value))
  heatmap_of_specificity_plot_2d <- heatmap_of_specificity_plot_2d + geom_tile(
    color = "white",
    lwd = 0.4,
    linetype = 1
  )
  heatmap_of_specificity_plot_2d <- heatmap_of_specificity_plot_2d + geom_text(
    aes(label = value),
    color = "white",
    size = 3
  )
  heatmap_of_specificity_plot_2d <- heatmap_of_specificity_plot_2d + labs(title = titlos)
  ## heatmapap_of_specificity_plot_2d <- heatmap_of_specificity_plot_2d + theme(
  ##   axis.title.x = element_blank(),
  ##   axis.title.y = element_blank()
  ## )
  heatmap_of_specificity_plot_2d <- heatmap_of_specificity_plot_2d + theme(
      axis.title = element_blank()
    )
  ## heatmap_of_specificity_plot_2d <- heatmap_of_specificity_plot_2d + theme_minimal()
  heatmap_of_specificity_plot_2d <- heatmap_of_specificity_plot_2d + scale_fill_gradient2( low = "#56B4E9", high = "#E69F00", midpoint = 0.5, limits = c(0.0, 1.0))
  heatmap_of_specificity_plot_2d <- heatmap_of_specificity_plot_2d + scale_y_discrete( drop = FALSE )
  heatmap_of_specificity_plot_2d <- heatmap_of_specificity_plot_2d + scale_x_discrete( drop = FALSE )

  return( heatmap_of_specificity_plot_2d )
  
}
