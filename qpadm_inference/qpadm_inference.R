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

# * Assignment of individuals to populations 

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

f2_blocks_for_single_model <- f2_from_geno(
  pref = input_prefixes[1],
  pops = individuals_in_populations,
  inds = individual_names
)
f2_blocks_for_single_model

## qpadm_multi iterates over a tibble whose rows represent
## right and left populations and target
all_qpadms <- qpadm_multi(
  data = f2_blocks_for_single_model,
  models = qp_models
)

lazaridis_battlegrounds <- function( result_of_qpadm_multi ){
  
  
}


### Following steps
## Make Lazaridis scheme DONE
## Run qpadm DONE
## Implement scoring function
## Loop previous 3
## Plot results
