# * Command line arguments

# * Assignment of individuals to populations

## Each of the 9 populations have 5 samples.
## Each of the outgroup populations have 10 samples.
## Outgroup samples are located after the ingroup samples in the vcf


normal_pops <- rep( sapply( 1:9, function(x){
  paste0( "pop_", x )
}), each = 5)

outgroups <- rep( sapply( 1:2, function(x){
  paste0( "outpop_", x )
}), each = 10)


all_names <- c(normal_pops, outgroups)
