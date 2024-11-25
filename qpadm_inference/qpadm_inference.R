# * Libraries

list_of_packages <- c(
  "ggplot2", "devtools",
  "argparse", "stringr",
  "Cairo", "tibble"
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

input_files <- list.files(path = input_folder, pattern = ".geno", full.names = TRUE)

## Debug ##
input_files <- list.files(
  path = "/media/storage/stef_sim/inference_estimation/sequencies/no_migration_constant_size/eig",
  pattern = ".geno", full.names = TRUE)
###########

## Need to take the full path of the prefixes to pass to f2_blocks
input_prefixes <- gsub(pattern = ".geno", replacement = "", input_files)

# * Data Preprocessing

# ** Assignment of individuals to populations 

## Each of the 9 populations have 5 samples.
## Each of the outgroup populations have 10 samples.
## Outgroup samples are located after the ingroup samples in the vcf

normal_pops <- rep( sapply( 0:8, function(x){
  paste0( "pop_", x )
}), each = 5)

outgroups <- rep(sapply(0:1, function(x) {
  paste0("outpop_", x)
}), each = 10)

## f2_blocks needs a vector with individual names and another vector
## of equal size that specifies the population of each individual.
individuals_in_populations <- c(normal_pops, outgroups)
population_names <- unique(individuals_in_populations)
individual_names <- paste("tsk_", 0:64, "indv", sep = '')

target <- "pop_4"
all_ancestors <- paste("pop_", c(0:3,5:8), sep = '')

# ** Build all models

## Following lines take all unique pairs of relatives.
## E.G. combination (pop_0, pop_1) exists but not (pop_1, pop_0)
all_models <- c()
for( i in 1:(length(all_ancestors)-1) ){
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
for (i in 1:nrow(all_models)) {
  right_pops <- population_names[!(population_names %in% c(all_models[i, ], target))]
  ## `list_of_right_pops`  will include all the permutations of right populations where only 1 is missing
  list_of_right_pops <- lapply(1:(length(right_pops) - 2), FUN = function(x) {
    right_pops[-x]
  })
  names(list_of_right_pops) <- right_pops[1:(length(right_pops) - 2)]
  list_of_right_pops$all <- right_pops  ## Lastly, add the model where none of the populations are missing
  for (j in 1:length(list_of_right_pops)) {
    left_all <- c(left_all, list(all_models[i, ]))
    right_all <- c(right_all, list(list_of_right_pops[[j]]))
    target_all <- c(target_all, target)
    exclude <- c(exclude, names(list_of_right_pops)[[j]])
  }
}

qp_models <- tibble(
  left = left_all,
  right = right_all,
  target = target_all
)

## Need to loop the following for all replicates

# * Analysis

# ** Run F2 and qpAdm

f2_blocks_for_single_model <- f2_from_geno(
  pref = input_prefixes[1],
  pops = individuals_in_populations,
  inds = individual_names
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

for( i in 1:length(all_qpadms) ){
    feasibility <- c(feasibility, all_qpadms[[i]]$popdrop$feasible[1])
    weights[[i]] <- unname(all_qpadms[[i]]$popdrop[1,c(7,8)])
    p_values <- c(p_values, all_qpadms[[i]]$popdrop[1,5])
}

## `data_v1` will be the data frame with all the information tidied up

data_v1 <- data.frame(
  left_1 = do.call(rbind, left_all)[,1],
  left_2 = do.call(rbind, left_all)[,2],
  weights_1 = as.data.frame(do.call( rbind, weights ))[,1],
  weights_2 = as.data.frame(do.call( rbind, weights ))[,2],
  exclude = exclude,
  p_values = unname(unlist(p_values)),
  p_values_all = unname(unlist(p_values[rep( which(exclude == 'all'), each = 7 )])),
  feasible = feasibility,
  feasible_all = feasibility[rep( which(exclude == 'all'), each = 7 )]
)

p_values_005 <- data_v1$p_values < 0.05
p_values_all_005 <- data_v1$p_values_all < 0.05

# ** Scoring loop

data_for_scoring_function <- as.data.frame( data_v1[exclude != "all",c(1,2,5,6,7,8,9)],
                                           row.names = 1:sum(exclude != "all") )

## Scoring function as Lazaridis described
scores <- c()
for( i in 1:nrow(data_for_scoring_function) ){
  if( data_for_scoring_function[i,]$p_values_all >= 0.05 ){
    if( data_for_scoring_function[i,]$p_values >= 0.05 ){
      ## Inference keeps being good regardless of inclusion.
      scores <- c(scores, 1)
      next
    } else {
      ## Inference is good when population was included but bad when it was excluded.
      scores <- c(scores, -1)
      next
    }
  } else {
    if( data_for_scoring_function[i,]$p_values >= 0.05 ){
      ## Inference was bad when population was included but good when it was excluded.
      ## This was not included in the original Supplementary.
      scores <- c(scores, -1)
      next
    } else {
      ## If inference is bad regardless of inclusion
      scores <- c(scores, 0)
      next
    }
  }
}

data_for_scoring_function$scores <- scores

## Scoring function as i think it is better
score_stef <- c()
update_score <- 0
for( i in 1:nrow(data_for_scoring_function) ){
  if( (data_for_scoring_function$p_values[i] < 0.05) & (data_for_scoring_function$p_values_all[i] < 0.05)  ){
    ## Both models are bad
    score_stef <- c( score_stef, 0 )
    next
  }else if( (data_for_scoring_function$p_values[i] < 0.05) & (data_for_scoring_function$p_values_all[i] > 0.05) ){
    ## Include is better, so I need to check if it is also feasible.
    if( data_for_scoring_function$feasible_all[i] == TRUE ){
      score_stef <- c( score_stef, 1 )
      next
    }else {
      score_stef <- c( score_stef, 0 )
      next
    }
  }else if( (data_for_scoring_function$p_values[i] > 0.05) & (data_for_scoring_function$p_values_all[i] < 0.05) ){
    ## Exclude is better, so I need to check if it is also feasible.
    if( data_for_scoring_function$feasible[i] == TRUE ){
      score_stef <- c( score_stef, -1 )
      next
    }else {
      score_stef <- c( score_stef, 0 )
      next
    }
  }else{ ## Both P_include and P_exclude are > 0.05. Need to check if they are also feasible.
    if( (data_for_scoring_function$feasible[i] == FALSE) & (data_for_scoring_function$feasible_all[i] == FALSE) ){
      ## Both are not feasible => tie.
      score_stef <- c(score_stef, 0)
      next
    }else if( (data_for_scoring_function$feasible[i] == FALSE) & (data_for_scoring_function$feasible_all[i] == TRUE) ){
      score_stef <- c(score_stef, 1) ## Need to discuss this. Include > Exclude
      next
    }else if( (data_for_scoring_function$feasible[i] == TRUE) & (data_for_scoring_function$feasible_all[i] == FALSE) ){
      score_stef <- c(score_stef, -1) # Exclude > Include
      next
    }else{ ## Both P-values are > 0.05 and both are feasible => Inference is good
      score_stef <- c(score_stef, 1)
      next
    }
  }
}
data_for_scoring_function$score_stef <- score_stef

### Following steps
## Make Lazaridis scheme DONE
## Run qpadm DONE
## Implement scoring function
## Loop previous 3
## Plot results
