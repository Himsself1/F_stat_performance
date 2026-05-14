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

if(!require("admixtools", character.only = TRUE)) {
  pak::pak("uqrmaie1/admixtools")
  ## devtools::install_github("uqrmaie1/admixtools")
}

read_table2 <- readr::read_table

# * Debug options

## input_folder <- "/media/storage/stef_sim/inference_estimation/harney_model_mut_1.25e-8_rec_1.25e-8_seq_249e+6/stats"
## plot_folder <- "/media/storage/stef_sim/inference_estimation/harney_model_mut_1.25e-8_rec_1.25e-8_seq_249e+6/plots"
## plot_functions <- "/home/stefanos/F_stat_performance/qpadm_inference/best_populations_plot_funtions.R"
## name_prefix <- "harney_model_mut_1.25e-8_rec_1.25e-8_seq_249e+6"


# * Command line arguments

args<-commandArgs(TRUE)
input_folder <- args[1]
plot_folder <- args[2]
plot_functions <- args[3]
name_prefix <- args[4]

print(input_folder)
print(plot_folder)
print(plot_functions)
print(name_prefix)

source(plot_functions)

dir.create(plot_folder, recursive = T)
input_names <- list.files(input_folder, pattern = ".tsv$", full.names = TRUE)

# * Data processing

all_ancestors <- factor(
  paste("pop_", c(0,3,5,6,7,8,9,10,12), sep = ''),
  levels = paste("pop_", c(0,3,5,6,7,8,9,10,12), sep = ''),
  ordered = T)

list_of_all_summaries <- list(
  accepted_models = list(),
  accepted_models_2d = list(),
  pop_scores = list(),
  pop_scores_2D = list(),
  specificity = list(),
  specificity_2d = list(),
  model_results = list(),
  good_models = list()
)

for (rep in 1:length(input_names)) {

  data_v1 = read.table(input_names[rep], header = T, sep = '\t')
  data_v1$left_1 <- factor(data_v1$left_1, levels = all_ancestors, ordered = T)
  data_v1$left_2 <- factor(data_v1$left_2, levels = all_ancestors, ordered = T)
  lapply(1:nrow(data_v1), function(x) {
    if(data_v1$left_1[x] > data_v1$left_2[x]) {
      data_v1[x, c(1,2)] <- data_v1[x, c(2,1)]
      }
    }) ## Check if order of populations is correct
  list_of_all_summaries$model_results[[rep]] <- data_v1

# ** Calculate "good" models per source.
  ## "good" are the models whose p_value is > 0.05 AND are "feasible"
  ## Sometimes qpAdm outputs weights not in [0,1]. These are not "feasible" models
  exclude <- data_v1$exclude
  factor_models <- list_of_all_summaries$model_results[[rep]][, 1:2] ## collumns of left_1 and left_2
  good_models <- list_of_all_summaries$model_results[[rep]][(list_of_all_summaries$model_results[[rep]]$feasible == TRUE) & (list_of_all_summaries$model_results[[rep]]$p_values > 0.05), ]
  list_of_all_summaries$good_models[[rep]] <- nrow(good_models)
  good_pops <- c(good_models$left_1, good_models$left_2)
  good_pops_2D <- data.frame(
    left_1 = good_models$left_1,
    left_2 = good_models$left_2
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
    filter(list_of_all_summaries$model_results[[rep]][,c(1, 2, 5, 6, 7, 8, 9)], exclude != "all"),
    row.names = 1:sum(exclude != "all")
  )
  
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
  results_versus <- data.frame()
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

  results_versus$left <- all_ancestors[results_versus$left]
  results_versus$exclude <- all_ancestors[results_versus$exclude]
  results_versus$direction <- as.numeric(results_versus$direction)
  ## I have upper triangular population pairs (if I have A & B, I don't calculate B & A).
  ## In order to score each population seperately, I need to concatenate "left"+"exclude" and
  ## flip the sign of the "exclude" population.
  results_versus_2nd_dim <- data.frame(
    left = c(results_versus[, 1], results_versus[, 2]),
    score = c(results_versus[, 3], results_versus[, 3] * (-1))
  )

  ## This shows how many times a population is the winner

  pop_scores <- as.data.frame(results_versus_2nd_dim %>%
                               group_by(left) %>%
                               summarise(score = sum(score)) %>%
                               arrange(desc(score)))

  pop_scores_2D <- as.data.frame(results_versus %>%
                                 ## group_by(left, exclude) %>%
                                 summarise(score = sum(direction), .by = c(left,exclude)))

  list_of_all_summaries$pop_scores[[rep]] <- pop_scores
  list_of_all_summaries$pop_scores_2D[[rep]] <- pop_scores_2D
}

# * Plot names

best_population_plot_name <- paste0(
  c( plot_folder, "/",
    name_prefix, "_best_pops.pdf"),
  collapse = "")

best_population_pair_plot_name <- paste0(
  c( plot_folder, "/",
    name_prefix, "_best_pop_pair.pdf"),
  collapse = "")

accepted_models_plot_name <- paste0(
  c( plot_folder, "/",
    name_prefix, "_accepted_models.pdf"),
  collapse = "")

accepted_models_2d_plot_name <- paste0(
  c( plot_folder, "/",
    name_prefix, "_accepted_models_2d.pdf"),
  collapse = "")

specificity_plot_name <- paste0(
  c( plot_folder, "/",
    name_prefix, "_specificity.pdf"),
  collapse = "")

specificity_2d_plot_name <- paste0(
  c( plot_folder, "/",
    name_prefix, "_specificity_2d.pdf"),
  collapse = "")

# * Calling functions & Printing plots

## Accepted Models

barplot_of_accepted_models <- plot_accepted_models(
  list_of_all_summaries$accepted_models,
  list_of_all_summaries$good_models,
  all_ancestors,
  ""
)

CairoPDF(accepted_models_plot_name)
barplot_of_accepted_models
dev.off()

heatmap_of_accepted_models_2d <- plot_accepted_models_2d(
  list_of_all_summaries$accepted_models_2d,
  list_of_all_summaries$good_models,
  all_ancestors,
  ""
)

CairoPDF( accepted_models_2d_plot_name )
heatmap_of_accepted_models_2d
dev.off()

## Best Populations

barplot_of_best_populations <- plot_best_2_pops(
  list_of_all_summaries$pop_scores,
  list_of_all_summaries$accepted_models,
  list_of_all_summaries$good_models,
  all_ancestors,
  ""
)

CairoPDF( best_population_plot_name )
barplot_of_best_populations
dev.off()

heatmap_of_best_pop_pair <- plot_best_pop_pair(
  list_of_all_summaries$pop_scores,
  list_of_all_summaries$accepted_models,
  list_of_all_summaries$good_models,
  all_ancestors,
  ""
)

CairoPDF( best_population_pair_plot_name )
heatmap_of_best_pop_pair
dev.off()

## Specificity

barplot_of_specificity <- plot_specificity(
  list_of_all_summaries$specificity,
  list_of_all_summaries$good_models,
  all_ancestors,
  ""
)

CairoPDF( specificity_plot_name )
barplot_of_specificity
dev.off()

heatmap_of_specificity_2d <- plot_specificity_2d(
  list_of_all_summaries$specificity_2d,
  list_of_all_summaries$good_models,
  all_ancestors,
  ""
)

CairoPDF( specificity_2d_plot_name )
heatmap_of_specificity_2d
dev.off()
