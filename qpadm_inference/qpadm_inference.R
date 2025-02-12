# * Libraries

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

# * Command line arguments

args<-commandArgs(TRUE)
input_folder <- args[1]
output_folder <- args[2]
plot_functions <- args[3]

## Debugging purposes ##
## input_folder <- "/media/storage/stef_sim/inference_estimation/sequencies/no_migration_constant_size/eig"
## output_folder <- "/media/storage/stef_sim/inference_estimation/statistics/no_migration_constant_size"

input_folder <- "/media/storage/stef_sim/inference_estimation/sequencies/migration_0123to5678_mig_02_constant_size_rec_e8/eig"
plot_functions <- "/home/stefanos/F_stat_performance/qpadm_inference/best_populations_plot_funtions.R"

########################

## Functions for plotting are located in another file
source( plot_functions )

input_files <- list.files(path = input_folder, pattern = ".geno", full.names = TRUE)
snp_files <- list.files(path = input_folder, pattern = ".snp", full.names = TRUE)

## parent_folder <- dirname(input_folder)

input_prefixes <- gsub(pattern = ".geno", replacement = "", input_files)
parent_dir <- dirname(input_folder)
metadata_file <- list.files(path = parent_dir, pattern = ".tsv", full.names = TRUE)

base_dir <- basename(parent_dir) ## Name of the parent folder of the simulation.
output_folder_for_plots <- paste0(output_folder, "/plots/", collapse = '')
output_folder_for_stats <- paste0(output_folder, "/stats/", collapse = '')
dir.create( output_folder_for_plots, recursive = TRUE )
dir.create( output_folder_for_stats, recursive = TRUE )

# * Data Preprocessing

# ** Assignment of individuals to populations 

metadata_info <- read.table( metadata_file, header = T, sep = '\t' )

## Each of the 9 populations have 5 samples.
## Each of the outgroup populations have 10 samples.
## Outgroup samples are located after the ingroup samples in the vcf

normal_pops <- rep( sapply( 0:8, function(x){
  paste0( "pop_", x )
}), each = 10)

outgroups <- rep(sapply(1:0, function(x) {
  paste0("outpop_", x)
}), each = 10)

## outgroups <- rep(sapply(0, function(x) {
##   paste0("outpop_", x)
## }), each = 10)


## f2_blocks needs a vector with individual names and another vector
## of equal size that specifies the population of each individual.
individuals_in_populations <- c(normal_pops, outgroups)
population_names <- unique(metadata_info$Population)
individual_names <- metadata_info$Ind_ID

target <- "pop_4"
all_ancestors <- paste("pop_", c(0:3,5:8), sep = '')

# ** Build all models

## Following lines take all unique pairs of relatives.
## E.G. combination (pop_0, pop_1) exists but not (pop_1, pop_0)
all_models <- c()
for( i in 1:(length(all_ancestors)-1) ) {
  for( j in (i+1):length(all_ancestors) ){
    all_models <- rbind( all_models, c(all_ancestors[i], all_ancestors[j]) )
  }
}

## Need to make a structure that contains all the combinations of `right` populations,
## as per Lazaridis (2024) scheme, i.e. dropping one `right` population each time.
## Combining each index of the variables bellow will consist of a qpadm model.
left_all <- list()
right_all <- list()
target_all <- list()
exclude <- c()

## Use `population_names` to control the right populations
for (i in 1:nrow(all_models)) {
  right_pops <- population_names[!(population_names %in% c(all_models[i, ], target))]
  ## `list_of_right_pops`  will include all the permutations of right populations where only 1 is missing
  ## We use -1 because there are 1 outgroups
  list_of_right_pops <- lapply(1:(length(right_pops) - length(unique(outgroups))), FUN = function(x) {
    right_pops[-x]
  })
  names(list_of_right_pops) <- right_pops[1:(length(right_pops) - length(unique(outgroups)))]
  list_of_right_pops$all <- right_pops  ## Lastly, add the model where none of the populations are missing
  for (j in 1:length(list_of_right_pops)) {
    left_all <- c(left_all, list(all_models[i, ]))
    right_all <- c(right_all, list(list_of_right_pops[[j]]))
    target_all <- c(target_all, target)
    exclude <- c(exclude, names(list_of_right_pops)[[j]])
  }
}

## I dont't need to be calculating qpadm models every iteration
qp_models <- tibble(
  left = left_all,
  right = right_all,
  target = target_all
)


# * Analysis

## This list will have all the info
list_of_all_summaries <- list(
  accepted_models = list(),
  accepted_models_2d = list(),
  best_pops = list(),
  best_pops_2D = list(),
  specificity = list(),
  specificity_2d = list()
)

for (rep in 1:length(input_prefixes)) {
# ** Run F2 and qpAdm
  
  f2_blocks_for_single_model <- f2_from_geno(
    pref = input_prefixes[rep],
    pops = individuals_in_populations,
    inds = individual_names,
    adjust_pseudohaploid = TRUE,
    blgsize = 0.05
  )
  
  ## qpadm_multi iterates over a tibble whose rows represent
  ## right and left populations and target
  all_qpadms <- qpadm_multi(
    data = f2_blocks_for_single_model,
    models = qp_models
  )
  
# ** Summarize the data

  feasibility <- c()
  weights <- list()
  p_values <- c()
  for (i in 1:length(all_qpadms)) {
    feasibility <- c(feasibility, all_qpadms[[i]]$popdrop$feasible[1])
    weights[[i]] <- unname(all_qpadms[[i]]$popdrop[1, c(7, 8)])
    p_values <- c(p_values, all_qpadms[[i]]$popdrop[1, 5])
  }
  
  ## `data_v1` will be the data frame with all the information tidied up
  ## In the final data frame, `p_value_all` and `feasible_all` show the p_value and
  ## feasibility of the model that doesn't exclude the respective population.
  data_v1 <- data.frame(
    left_1 = factor(do.call(rbind, left_all)[, 1], levels = all_ancestors),
    left_2 = factor(do.call(rbind, left_all)[, 2], levels = all_ancestors),
    weights_1 = as.data.frame(do.call(rbind, weights))[, 1],
    weights_2 = as.data.frame(do.call(rbind, weights))[, 2],
    exclude = exclude,
    p_values = unname(unlist(p_values)),
    p_values_all = unname(unlist(p_values[rep(which(exclude == "all"), each = length(all_ancestors) - 1)])),
    feasible = feasibility,
    feasible_all = feasibility[rep(which(exclude == "all"), each = length(all_ancestors) - 1)],
    replicate = rep( rep, length(exclude) )
  )
  
# ** Calculate "good" models per source. 
  ## "good" are the models whose p_value is > 0.05 AND are "feasible"
  ## Sometimes qpAdm outputs weights not in [0,1]. These are not "feasible" models
  
  factor_models <- data_v1[, 1:2] ## collumns of left_1 and left_2
  good_models <- data_v1[(data_v1$feasible == TRUE) & (data_v1$p_values > 0.05), ]
  good_pops <- factor(c(good_models$left_1, good_models$left_2), levels = all_ancestors)
  good_pops_2D <- data_frame(
    left_1 = factor(good_models$left_1, levels = all_ancestors),
    left_2 = factor(good_models$left_2, levels = all_ancestors)
  )
  
  percent_of_accepted_models <- table(good_pops) / choose(length(unique(exclude)), 2)
  percent_of_accepted_models_2d <- table(good_pops_2D) / table(factor_models)
  
  specificity <- table(good_pops) / length(good_pops)
  specificity_2d <- table(good_pops_2D) / length(good_pops_2D)
    
  list_of_all_summaries$accepted_models[[rep]] <- percent_of_accepted_models
  list_of_all_summaries$accepted_models_2d[[rep]] <- percent_of_accepted_models_2d

  list_of_all_summaries$specificity[[rep]] <- specificity
  list_of_all_summaries$specificity_2d[[rep]] <- specificity_2d
  
# ** Scoring loop
  
  data_for_scoring_function <- as.data.frame(
    data_v1[exclude != "all", c(1, 2, 5, 6, 7, 8, 9)],
    row.names = 1:sum(exclude != "all")
  )
  
  ## scores <- c()
  ## for( i in 1:nrow(data_for_scoring_function) ){
  ##   if( data_for_scoring_function[i,]$p_values_all >= 0.05 ){
  ##     if( data_for_scoring_function[i,]$p_values >= 0.05 ){
  ##       ## Inference keeps being good regardless of inclusion.
  ##       scores <- c(scores, 1)
  ##       next
  ##     } else {
  ##       ## Inference is good when population was included but bad when it was excluded.
  ##       scores <- c(scores, -1)
  ##       next
  ##     }
  ##   } else {
  ##     if( data_for_scoring_function[i,]$p_values >= 0.05 ){
  ##       ## Inference was bad when population was included but good when it was excluded.
  ##       ## This was not included in the original Supplementary.
  ##       scores <- c(scores, -1)
  ##       next
  ##     } else {
  ##       ## If inference is bad regardless of inclusion
  ##       scores <- c(scores, 0)
  ##       next
  ##     }
  ##   }
  ## }
  ## data_for_scoring_function$scores <- scores

# **  Calculate pairwise "resillience" (as described by Lazaridis 2024).
  ## Returns 1 if resillient, -1 if not resillient and 0 if a decision cannot be made.
  score_stef <- c()
  for (i in 1:nrow(data_for_scoring_function)) {
    if ((data_for_scoring_function$p_values[i] < 0.05) & (data_for_scoring_function$p_values_all[i] < 0.05)) {
      ## Both models are bad
      score_stef <- c(score_stef, 0)
      next
    } else if ((data_for_scoring_function$p_values[i] < 0.05) & (data_for_scoring_function$p_values_all[i] > 0.05)) {
      ## Include is better, so I need to check if it is also feasible.
      if (data_for_scoring_function$feasible_all[i] == TRUE) {
        score_stef <- c(score_stef, 0)
        next
      } else {
        score_stef <- c(score_stef, 0)
        next
      }
    } else if ((data_for_scoring_function$p_values[i] > 0.05) & (data_for_scoring_function$p_values_all[i] < 0.05)) {
      ## Exclude is better, so I need to check if it is also feasible.
      if (data_for_scoring_function$feasible[i] == TRUE) {
        score_stef <- c(score_stef, -1)
        next
      } else {
        score_stef <- c(score_stef, 0)
        next
      }
    } else {
      ## Both P_include and P_exclude are > 0.05. Need to check if they are also feasible.
      if ((data_for_scoring_function$feasible[i] == FALSE) & (data_for_scoring_function$feasible_all[i] == FALSE)) {
        ## Both are not feasible => tie.
        score_stef <- c(score_stef, 0)
        next
      } else if ((data_for_scoring_function$feasible[i] == FALSE) & (data_for_scoring_function$feasible_all[i] == TRUE)) {
        score_stef <- c(score_stef, 0) # Need to discuss this. Include > Exclude
        next
      } else if ((data_for_scoring_function$feasible[i] == TRUE) & (data_for_scoring_function$feasible_all[i] == FALSE)) {
        score_stef <- c(score_stef, -1) # Exclude > Include
        next
      } else {
        ## Both P-values are > 0.05 and both are feasible => Inference is good
        score_stef <- c(score_stef, 1)
        next
      }
    }
  }
  data_for_scoring_function$score_stef <- score_stef
  
  data_for_single_sim_heatmap <- data_for_scoring_function[, c(1:3, 8)]
  data_heatmap_2 <- data_for_scoring_function[, c(2, 1, 3, 8)]
  names(data_heatmap_2) <- names(data_for_single_sim_heatmap)
  data_versus <- rbind(data_for_single_sim_heatmap, data_heatmap_2)

# ** Calculating Lazaridis score

  ## If A is resillient to B and B is not resillient to A then A vs B = 1
  ## If A is resillient to B and B is resillient to A then A vs B = 0
  ## If A is not resillient to B and B is resillient to A then A vs B = -1
  ## The sign of the score denotes which population is favored in the comparison.
  results_versus <- data_frame()
  for (l1 in 1:(length(all_ancestors) - 1)) {
    for (ex in (l1 + 1):(length(all_ancestors))) {
      if (l1 == ex) {
        next
      }
      for (l2 in 1:length(all_ancestors)) {
        if (l2 %in% c(l1, ex)) {
          next
        } else {
          index_model_A <- which(
          (data_versus$left_1 == all_ancestors[l1]) &
            (data_versus$exclude == all_ancestors[ex]) &
            (data_versus$left_2 == all_ancestors[l2])
          )
          index_model_B <- which(
          (data_versus$left_1 == all_ancestors[ex]) &
            (data_versus$exclude == all_ancestors[l1]) &
            (data_versus$left_2 == all_ancestors[l2])
          )
          if (data_versus$score_stef[index_model_A] > data_versus$score_stef[index_model_B]) {
            temp_result_versus <- c(all_ancestors[l1], all_ancestors[ex], 1)
          } else if (data_versus$score_stef[index_model_A] < data_versus$score_stef[index_model_B]) {
            temp_result_versus <- c(all_ancestors[l1], all_ancestors[ex], -1)
          } else {
            temp_result_versus <- c(all_ancestors[l1], all_ancestors[ex], 0)
          }
        }
        results_versus <- rbind(results_versus, temp_result_versus)
      }
    }
  }
  names(results_versus) <- c("left", "exclude", "direction")

  results_versus$left <- factor(results_versus$left, levels = all_ancestors)
  results_versus$exclude <- factor(results_versus$exclude, levels = all_ancestors)
  results_versus$direction <- as.numeric(results_versus$direction)
  ## I have upper triangular population pairs (if I have A & B, I don't calculate B & A).
  ## In order to score each population seperately, I need to concatenate "left"+"exclude" and
  ## flip the sign of the "exclude" population.
  results_versus_2nd_dim <- data_frame(
    left = factor(c(results_versus[, 1], results_versus[, 2]), levels = all_ancestors),
    score = c(results_versus[, 3], results_versus[, 3] * (-1))
  )

  ## This shows how many times a population is the winner

  best_pops <- as.data.frame(results_versus_2nd_dim %>%
                               group_by(left) %>%
                               summarise(score = sum(score)) %>%
                               arrange(desc(score)))

  best_pops_2D <- as.data.frame(results_versus %>%
                                  group_by(left, exclude) %>%
                                  summarise(score = as.numeric(sum(direction))))

  list_of_all_summaries$best_pops[[rep]] <- best_pops
  list_of_all_summaries$best_pops_2D[[rep]] <- best_pops_2D
}

# * Plotting

# ** Names of files

best_population_plot_name <- paste0(
  c( output_folder_for_plots, "/",
    base_dir, "_best_pops.pdf"),
  collapse = "")

best_population_pair_plot_name <- paste0(
  c( output_folder_for_plots, "/",
    base_dir, "_best_pop_pair.pdf"),
  collapse = "")

accepted_models_plot_name <- paste0(
  c( output_folder_for_plots, "/",
    base_dir, "_accepted_models.pdf"),
  collapse = "")

accepted_models_2d_plot_name <- paste0(
  c( output_folder_for_plots, "/",
    base_dir, "_accepted_models_2d.pdf"),
  collapse = "")

specificity_plot_name <- paste0(
  c( output_folder_for_plots, "/",
    base_dir, "_specificity.pdf"),
  collapse = "")

specificity_2d_plot_name <- paste0(
  c( output_folder_for_plots, "/",
    base_dir, "_specificity_2d.pdf"),
  collapse = "")

# ** Calling functions & Printing plots

barplot_of_accepted_models <- plot_accepted_models(
  list_of_all_summaries$accepted_models,
  all_ancestors,
  ""
)

CairoPDF( accepted_models_plot_name )
barplot_of_accepted_models
dev.off()

heatmap_of_accepted_models_2d <- plot_accepted_models_2d(
  list_of_all_summaries$accepted_models_2d,
  all_ancestors,
  ""
)

## CairoPDF( "test_accepted.pdf" )
## heatmap_of_accepted_models_2d
## dev.off()

CairoPDF( accepted_models_2d_plot_name )
heatmap_of_accepted_models_2d
dev.off()

barplot_of_best_populations <- plot_best_2_pops(
  list_of_all_summaries$best_pops,
  all_ancestors,
  ""
)

CairoPDF( best_population_plot_name )
barplot_of_best_populations
dev.off()

heatmap_of_best_pop_pair <- plot_best_pop_pair(
  list_of_all_summaries$best_pops,
  all_ancestors,
  ""
)

CairoPDF( best_population_pair_plot_name )
heatmap_of_best_pop_pair
dev.off()

barplot_of_specificity <- plot_specificity(
  list_of_all_summaries$specificity,
  all_ancestors,
  ""
)

## CairoPDF( "specificity_plot.pdf" )
## barplot_of_specificity
## dev.off()

CairoPDF( specificity_plot_name )
barplot_of_specificity
dev.off()

heatmap_of_specificity_2d <- plot_specificity_2d(
  list_of_all_summaries$specificity_2d,
  all_ancestors,
  ""
)

## CairoPDF( "specificity_plot_2d.pdf" )
## heatmap_of_specificity_2d
## ## heatmap_of_specificity_2d
## dev.off()

CairoPDF( specificity_2d_plot_name )
heatmap_of_specificity_2d
dev.off()

## output_folder <- "/media/storage/stef_sim/inference_estimation/plots"

# * Legacy (non-automated code)

## # ** Scale 1

## # *** 2 Best Pops

## ## % of simulation that "population" was among the 2 best.

## best_2_pops_barplot_name <- paste0(c( output_folder, "/best_2_pops.pdf" ), collapse = '')

## best_2_populations <- melt(do.call(rbind, lapply(1:rep, function(x) {
##         as.character(list_of_all_summaries$best_pops[[x]][1:2, 1])
## })))

## best_2_for_barplot <- data.frame(
##         onomata = factor(names(table(factor(best_2_populations[, 3], levels = all_ancestors))), levels = all_ancestors),
##         values = as.vector(table(best_2_populations[, 3]))
## )

## barplot_of_best_2_pops <- ggplot(best_2_for_barplot, aes(x = onomata, y = as.numeric(values)))
## barplot_of_best_2_pops <- barplot_of_best_2_pops + geom_bar(stat = "identity", fill = "steelblue")
## barplot_of_best_2_pops <- barplot_of_best_2_pops + geom_text(
##   aes(label = values),
##   vjust = -0.3, size = 4
## )
## barplot_of_best_2_pops <- barplot_of_best_2_pops + theme(
##         axis.title = element_blank()
## )
## barplot_of_best_2_pops <- barplot_of_best_2_pops + labs(
##         title = "Scale 0.5"
## )


## CairoPDF( best_2_pops_barplot_name )
## barplot_of_best_2_pops
## dev.off()

## # *** 2 Best Pops 2D

## ## % of simulation that "population pair" was among the best pair.

## heatmap_of_best_two_pops_2d_plot_name <- paste0(c( output_folder, "/best_2_pops_2d.pdf" ), collapse = '')

## best_two_pops_2d <- as.data.frame(do.call(rbind, lapply(1:rep, function(x) {
##   sort(as.character(list_of_all_summaries$best_pops[[x]][1:2, 1]))
## })))
## names(best_two_pops_2d) <- c("best_1", "best_2")

## best_two_pops_2d$best_1 <- factor(best_two_pops_2d$best_1, levels = all_ancestors)
## best_two_pops_2d$best_2 <- factor(best_two_pops_2d$best_2, levels = all_ancestors)

## melted_best_two_pops_2d <- melt(table(best_two_pops_2d) / 100)
## triangle <- as.character(melted_best_two_pops_2d$best_1) < as.character(melted_best_two_pops_2d$best_2)

## heatmap_of_best_two_pops_2d <- ggplot(melted_best_two_pops_2d[triangle, ], aes(x = best_1, y = best_2, fill = value))
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + geom_tile(
##         color = "white",
##         lwd = 0.4,
##         linetype = 1
## )
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + geom_text(
##   aes(label = value),
##   color = "white",
##   size = 3
## )
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + labs(title = "Scale 0.5")
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + theme(
##   axis.title.x = element_blank(),
##   axis.title.y = element_blank()
## )
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + scale_fill_gradient( low = "#56B4E9", high = "#E69F00")

## CairoPDF( heatmap_of_best_two_pops_2d_plot_name )
## heatmap_of_best_two_pops_2d
## dev.off()

## # *** Acceplted Models

## ## Number of models of "population" that are accepted divided by
## ## all the possible models that the population can take part in.

## average_accepted_models_plot_name <- paste0(
##   c(output_folder, "/average_accepted_models.pdf"),
##         collapse = ""
## )

## average_accepted_models <- colMeans( do.call(rbind, list_of_all_summaries$accepted_models) )
## average_accepted_models_for_barplot <- data.frame(
##         onomata = factor(names(average_accepted_models), levels = all_ancestors),
##         values = as.numeric(sprintf("%.4f", average_accepted_models))
## )

## barplot_of_accepted_models <- ggplot(average_accepted_models_for_barplot, aes(x = onomata, y = as.numeric(values)))
## barplot_of_accepted_models <- barplot_of_accepted_models + geom_bar(stat = "identity", fill = "steelblue")
## barplot_of_accepted_models <- barplot_of_accepted_models + geom_text(
##   aes(label = values),
##   vjust = -0.3, size = 4
## )
## barplot_of_accepted_models <- barplot_of_accepted_models + labs(
##         title = "Scale 0.5"
## )
## barplot_of_accepted_models <- barplot_of_accepted_models + theme(
##         axis.title = element_blank()
## )

## CairoPDF( average_accepted_models_plot_name )
## barplot_of_accepted_models
## dev.off()

## # *** Accepted Models 2D

## ## Number of models of "population pair" that are accepted divided by
## ## all the possible models that the population pair can take part in.

## heatmap_of_accepted_models_2D_plot_name <- paste0(
##         c(output_folder, "/accepted_models_2d.pdf"),
##         collapse = ""
## )

## melted_accepted_models_2d <- melt(Reduce("+", list_of_all_summaries$accepted_models_2d) / 100, na.rm = TRUE)
## melted_accepted_models_2d <- melted_accepted_models_2d[!is.nan(melted_accepted_models_2d$value), ]
## melted_accepted_models_2d$value  <- as.numeric(sprintf("%.4f", melted_accepted_models_2d$value))

## heatmap_of_accepted_models_plot_2D <- ggplot(melted_accepted_models_2d, aes( x = left_1, y = left_2, fill = value))
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + geom_tile(
##         color = "white",
##         lwd = 0.4,
##         linetype = 1
## )
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + geom_text(
##   aes(label = value),
##   color = "white",
##   size = 3
## )
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + labs(title = "Scale 0.5")
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + theme(
##   axis.title.x = element_blank(),
##   axis.title.y = element_blank()
## )
## ## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + theme_minimal()
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + scale_fill_gradient( low = "#56B4E9", high = "#E69F00")

## CairoPDF( heatmap_of_accepted_models_2D_plot_name )
## heatmap_of_accepted_models_plot_2D
## dev.off()

## # ** Scale 2

## # *** 2 Best Pops

## best_2_pops_barplot_name <- paste0(c( output_folder, "/best_2_pops_scale_2.pdf" ), collapse = '')

## best_2_populations <- melt(do.call(rbind, lapply(1:rep, function(x) {
##         as.character(list_of_all_summaries$best_pops[[x]][1:2, 1])
## })))

## best_2_for_barplot <- data.frame(
##   onomata = factor(names(table(factor(best_2_populations[, 3], levels = all_ancestors))), levels = all_ancestors),
##         values = as.vector(table(factor(best_2_populations[, 3], levels = all_ancestors)))
## )

## barplot_of_best_2_pops <- ggplot(best_2_for_barplot, aes(x = onomata, y = as.numeric(values)))
## barplot_of_best_2_pops <- barplot_of_best_2_pops + geom_bar(stat = "identity", fill = "steelblue")
## barplot_of_best_2_pops <- barplot_of_best_2_pops + geom_text(
##   aes(label = values),
##   vjust = -0.3, size = 4
## )
## barplot_of_best_2_pops <- barplot_of_best_2_pops + theme(
##         axis.title = element_blank()
## )
## barplot_of_best_2_pops <- barplot_of_best_2_pops + labs(
##         title = "Scale 1"
## )


## CairoPDF( best_2_pops_barplot_name )
## barplot_of_best_2_pops
## dev.off()

## # *** 2 Best Pops 2D

## heatmap_of_best_two_pops_2d_plot_name <- paste0(c( output_folder, "/best_2_pops_2d_scale_2.pdf" ), collapse = '')

## best_two_pops_2d <- as.data.frame(do.call(rbind, lapply(1:rep, function(x) {
##   sort(as.character(list_of_all_summaries$best_pops[[x]][1:2, 1]))
## })))
## names(best_two_pops_2d) <- c("best_1", "best_2")

## best_two_pops_2d$best_1 <- factor(best_two_pops_2d$best_1, levels = all_ancestors)
## best_two_pops_2d$best_2 <- factor(best_two_pops_2d$best_2, levels = all_ancestors)

## melted_best_two_pops_2d <- melt(table(best_two_pops_2d) / 100)
## triangle <- as.character(melted_best_two_pops_2d$best_1) < as.character(melted_best_two_pops_2d$best_2)

## heatmap_of_best_two_pops_2d <- ggplot(melted_best_two_pops_2d[triangle, ], aes(x = best_1, y = best_2, fill = value))
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + geom_tile(
##         color = "white",
##         lwd = 0.4,
##         linetype = 1
## )
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + geom_text(
##   aes(label = value),
##   color = "white",
##   size = 3
## )
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + labs(title = "Scale 1")
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + theme(
##   axis.title.x = element_blank(),
##   axis.title.y = element_blank()
## )
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + scale_fill_gradient( low = "#56B4E9", high = "#E69F00")

## CairoPDF( heatmap_of_best_two_pops_2d_plot_name )
## heatmap_of_best_two_pops_2d
## dev.off()

## # *** Acceplted Models

## average_accepted_models_plot_name <- paste0(
##   c(output_folder, "/average_accepted_models_scale_2.pdf"),
##         collapse = ""
## )

## average_accepted_models <- colMeans( do.call(rbind, list_of_all_summaries$accepted_models) )
## average_accepted_models_for_barplot <- data.frame(
##         onomata = factor(names(average_accepted_models), levels = all_ancestors),
##         values = as.numeric(sprintf("%.4f", average_accepted_models))
## )

## barplot_of_accepted_models <- ggplot(average_accepted_models_for_barplot, aes(x = onomata, y = as.numeric(values)))
## barplot_of_accepted_models <- barplot_of_accepted_models + geom_bar(stat = "identity", fill = "steelblue")
## barplot_of_accepted_models <- barplot_of_accepted_models + geom_text(
##   aes(label = values),
##   vjust = -0.3, size = 4
## )
## barplot_of_accepted_models <- barplot_of_accepted_models + labs(
##         title = "Scale 1"
## )
## barplot_of_accepted_models <- barplot_of_accepted_models + theme(
##         axis.title = element_blank()
## )

## CairoPDF( average_accepted_models_plot_name )
## barplot_of_accepted_models
## dev.off()

## # *** Accepted Models 2D

## heatmap_of_accepted_models_2D_plot_name <- paste0(
##         c(output_folder, "/accepted_models_2d_scale_2.pdf"),
##         collapse = ""
## )

## melted_accepted_models_2d <- melt(Reduce("+", list_of_all_summaries$accepted_models_2d) / 100, na.rm = TRUE)
## melted_accepted_models_2d <- melted_accepted_models_2d[!is.nan(melted_accepted_models_2d$value), ]
## melted_accepted_models_2d$value  <- as.numeric(sprintf("%.4f", melted_accepted_models_2d$value))

## heatmap_of_accepted_models_plot_2D <- ggplot(melted_accepted_models_2d, aes( x = left_1, y = left_2, fill = value))
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + geom_tile(
##         color = "white",
##         lwd = 0.4,
##         linetype = 1
## )
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + geom_text(
##   aes(label = value),
##   color = "white",
##   size = 3
## )
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + labs(title = "Scale 1")
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + theme(
##   axis.title.x = element_blank(),
##   axis.title.y = element_blank()
## )
## ## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + theme_minimal()
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + scale_fill_gradient( low = "#56B4E9", high = "#E69F00")

## CairoPDF( heatmap_of_accepted_models_2D_plot_name )
## heatmap_of_accepted_models_plot_2D
## dev.off()

## # ** Scale 5

## # *** 2 Best Pops

## best_2_pops_barplot_name <- paste0(c( output_folder, "/best_2_pops_scale_5.pdf" ), collapse = '')

## best_2_populations <- melt(do.call(rbind, lapply(1:rep, function(x) {
##         as.character(list_of_all_summaries$best_pops[[x]][1:2, 1])
## })))

## best_2_for_barplot <- data.frame(
##   onomata = factor(names(table(factor(best_2_populations[, 3], levels = all_ancestors))), levels = all_ancestors),
##         values = as.vector(table(best_2_populations[, 3]))
## )

## barplot_of_best_2_pops <- ggplot(best_2_for_barplot, aes(x = onomata, y = as.numeric(values)))
## barplot_of_best_2_pops <- barplot_of_best_2_pops + geom_bar(stat = "identity", fill = "steelblue")
## barplot_of_best_2_pops <- barplot_of_best_2_pops + geom_text(
##   aes(label = values),
##   vjust = -0.3, size = 4
## )
## barplot_of_best_2_pops <- barplot_of_best_2_pops + theme(
##         axis.title = element_blank()
## )
## barplot_of_best_2_pops <- barplot_of_best_2_pops + labs(
##         title = "Scale 2.5"
## )


## CairoPDF( best_2_pops_barplot_name )
## barplot_of_best_2_pops
## dev.off()

## # *** 2 Best Pops 2D


## heatmap_of_best_two_pops_2d_plot_name <- paste0(c( output_folder, "/best_2_pops_2d_scale_5.pdf" ), collapse = '')

## best_two_pops_2d <- as.data.frame(do.call(rbind, lapply(1:rep, function(x) {
##   sort(as.character(list_of_all_summaries$best_pops[[x]][1:2, 1]))
## })))
## names(best_two_pops_2d) <- c("best_1", "best_2")

## best_two_pops_2d$best_1 <- factor(best_two_pops_2d$best_1, levels = all_ancestors)
## best_two_pops_2d$best_2 <- factor(best_two_pops_2d$best_2, levels = all_ancestors)

## melted_best_two_pops_2d <- melt(table(best_two_pops_2d) / 100)
## triangle <- as.character(melted_best_two_pops_2d$best_1) < as.character(melted_best_two_pops_2d$best_2)

## heatmap_of_best_two_pops_2d <- ggplot(melted_best_two_pops_2d[triangle, ], aes(x = best_1, y = best_2, fill = value))
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + geom_tile(
##         color = "white",
##         lwd = 0.4,
##         linetype = 1
## )
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + geom_text(
##   aes(label = value),
##   color = "white",
##   size = 3
## )
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + labs(title = "Scale 2.5")
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + theme(
##   axis.title.x = element_blank(),
##   axis.title.y = element_blank()
## )
## heatmap_of_best_two_pops_2d <- heatmap_of_best_two_pops_2d + scale_fill_gradient( low = "#56B4E9", high = "#E69F00")

## CairoPDF( heatmap_of_best_two_pops_2d_plot_name )
## heatmap_of_best_two_pops_2d
## dev.off()

## # *** Acceplted Models

## average_accepted_models_plot_name <- paste0(
##   c(output_folder, "/average_accepted_models_scale_5.pdf"),
##         collapse = ""
## )

## average_accepted_models <- colMeans( do.call(rbind, list_of_all_summaries$accepted_models) )
## average_accepted_models_for_barplot <- data.frame(
##         onomata = factor(names(average_accepted_models), levels = all_ancestors),
##         values = as.numeric(sprintf("%.4f", average_accepted_models))
## )

## barplot_of_accepted_models <- ggplot(average_accepted_models_for_barplot, aes(x = onomata, y = as.numeric(values)))
## barplot_of_accepted_models <- barplot_of_accepted_models + geom_bar(stat = "identity", fill = "steelblue")
## barplot_of_accepted_models <- barplot_of_accepted_models + geom_text(
##   aes(label = values),
##   vjust = -0.3, size = 4
## )
## barplot_of_accepted_models <- barplot_of_accepted_models + labs(
##         title = "Scale 2.5"
## )
## barplot_of_accepted_models <- barplot_of_accepted_models + theme(
##         axis.title = element_blank()
## )

## CairoPDF( average_accepted_models_plot_name )
## barplot_of_accepted_models
## dev.off()

## # *** Accepted Models 2D

## heatmap_of_accepted_models_2D_plot_name <- paste0(
##         c(output_folder, "/accepted_models_2d_scale_5.pdf"),
##         collapse = ""
## )

## melted_accepted_models_2d <- melt(Reduce("+", list_of_all_summaries$accepted_models_2d) / 100, na.rm = TRUE)
## melted_accepted_models_2d <- melted_accepted_models_2d[!is.nan(melted_accepted_models_2d$value), ]
## melted_accepted_models_2d$value  <- as.numeric(sprintf("%.4f", melted_accepted_models_2d$value))

## heatmap_of_accepted_models_plot_2D <- ggplot(melted_accepted_models_2d, aes( x = left_1, y = left_2, fill = value))
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + geom_tile(
##         color = "white",
##         lwd = 0.4,
##         linetype = 1
## )
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + geom_text(
##   aes(label = value),
##   color = "white",
##   size = 3
## )
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + labs(title = "Scale 2.5")
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + theme(
##   axis.title.x = element_blank(),
##   axis.title.y = element_blank()
## )
## ## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + theme_minimal()
## heatmap_of_accepted_models_plot_2D <- heatmap_of_accepted_models_plot_2D + scale_fill_gradient( low = "#56B4E9", high = "#E69F00")

## CairoPDF( heatmap_of_accepted_models_2D_plot_name )
## heatmap_of_accepted_models_plot_2D
## dev.off()
