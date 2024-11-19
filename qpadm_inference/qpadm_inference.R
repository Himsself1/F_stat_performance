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

individuals_in_populations <- c(normal_pops, outgroups)
population_names <- unique(individuals_in_populations)

target <- "pop_4"
ancestors_of_pop3 <- paste("pop_", c(0:3), sep = '')
ancestors_of_pop5 <- paste("pop_", c(5:8), sep = '')
all_models <- expand.grid(ancestors_of_pop3, ancestors_of_pop5, stringsAsFactors = F)

## Runs f2 for all populations \
## TODO: Make a loop for all pairs of ancestors.
## Should I make pairs of ancestors from the same side of the tree?

f2_blocks_for_single_model <- f2_from_geno(
  pref = input_prefixes[1],
  inds = individuals_in_populations
)

## Need to make a structure that contains all the combinations of `right` populations.

## Each index of the variables bellow will consist of a qpadm model.
left_all <- list()
right_all <- list()
target_all <- list()

for (i in 1:nrow(all_models)) {
  right_pops <- population_names[!(population_names %in% c(all_models[i, ], target))]
  ## `list_of_right_pops ` will include all the permutations of right populations where only 1 is missing
  list_of_right_pops <- lapply(1:(length(right_pops) - 2), FUN = function(x) {
    right_pops[-x]
  })
  names(list_of_right_pops) <- right_pops[1:(length(right_pops) - 2)]
  list_of_right_pops$all <- right_pops  ## Lastly, add the model where none of the populations are missing
  for (j in 1:length(list_of_right_pops)) {
    left_all <- c(left_all, list(all_models[i, ]))
    right_all <- c(right_all, list(list_of_right_pops[[j]]))
    target_all <- c(target_all, target)
  }
}

qp_models <- tibble(
  left = left_all,
  right = right_all,
  target = target_all
)

?tibble

## right_pops <- population_names[!(population_names %in% all_models[1, ])]
## list_of_right_pops <- lapply(1:(length(right_pops) - 2), FUN = function(x) {
##   right_pops[-x]
## })
## names(list_of_right_pops) <- right_pops[1:(length(right_pops)-2)]
## list_of_right_pops$all <- right_pops

qpadm(
  f2_blocks_for_single_model,
  left = all_models[i, c(1, 3)],
  right = right_pops[[i]]
)

## Need to make a loop for qpadm for all sets of right pops given the ancestors

### Following steps
## Make Lazaridis scheme
## Run qpadm
## Implement scoring function
## Plot results
