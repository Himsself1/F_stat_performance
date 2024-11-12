# * Libraries

import msprime
import argparse
import pathlib

# * Command line arguments

cli = argparse.ArgumentParser(
    description = 'msprime simulations with ',
    prog = 'simulations_msprime.py'
)

cli.add_argument( "-out_folder", type = str, required = True,
                  help = "Full path to output folder. The program will create all child folders that may not exist.")
cli.add_argument( "-name", type = str, required = True,
                  help = "Name prefix of output vcf files.")
cli.add_argument( "-how_many", type = int, default = 1,
                  help = 'How many simulations are going to be run in one go.' )
cli.add_argument( "-mig_ratio", type = float, required = True,
                  help = 'Persentage of migrant in migration event.')
argue = cli.parse_args()

pathlib.Path(argue.out_folder).mkdir(parents = True, exist_ok = True) # Makes Output directories

# * Build population model 

# ** Configure split times

## These times are in generations

split_time_2_3 = 80
split_time_1_2_3 = 150
split_time_0_1_2_3 = 200
split_time_5_6 = 80
split_time_7_5_6 = 150
split_time_8_7_5_6 = 200
split_time_ingroups = 400
split_time_ingroups_out_2 = 500
split_time_ancestral_all = 1000

admixture_time = 40

# ** Configure Population Sizes

# pop_size_0 = 1000
# pop_size_1 = 1000
# pop_size_2 = 1000
# pop_size_3 = 1000
# pop_size_4 = 1000
# pop_size_5 = 1000
# pop_size_6 = 1000
# pop_size_7 = 1000
# pop_size_8 = 1000

pop_initial_size = [1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000]

# pop_size_out_1 = 1000
# pop_size_out_0 = 1000
outgroup_pop_size = [1000, 1000]

pop_size_ancestral_2_3 = 1000
pop_size_ancestral_1_2_3 = 1000
pop_size_ancestral_0_1_2_3 = 1000
pop_size_ancestral_5_6 = 1000
pop_size_ancestral_7_5_6 = 1000
pop_size_ancestral_8_7_5_6 = 1000
pop_size_ancestral_ingroups = 1000
pop_size_ancestral_ingroups_out_1 = 1000
pop_size_ancestral_all = 1000

demography = msprime.Demography()

# ** Create populations 0-9

[demography.add_population(name= "pop_"+str(i), initial_size = pop_initial_size[i] ) for i in range(9)]

# ** Create 2 outgroups for reference

[demography.add_population(name= "outpop_"+str(i), initial_size = outgroup_pop_size[i] ) for i in range(2)]

# ** Create ancestral populations

demography.add_population(name = "ancestral_2_3", initial_size = pop_size_ancestral_2_3)
demography.add_population(name = "ancestral_1_2_3", initial_size = pop_size_ancestral_1_2_3)
demography.add_population(name = "ancestral_0_1_2_3", initial_size = pop_size_ancestral_0_1_2_3)
demography.add_population(name = "ancestral_5_6", initial_size = pop_size_ancestral_5_6)
demography.add_population(name = "ancestral_7_5_6", initial_size = pop_size_ancestral_7_5_6)
demography.add_population(name = "ancestral_8_7_5_6", initial_size = pop_size_ancestral_8_7_5_6)
demography.add_population(name = "ancestral_ingroups", initial_size = pop_size_ancestral_ingroups)
demography.add_population(name = "ancestral_ingroups_out_1", initial_size = pop_size_ancestral_ingroups_out_1)
demography.add_population(name = "ancestral_all", initial_size = pop_size_ancestral_all)

# ** Create coalescent events

demography.add_population_split(
    time = split_time_2_3,
    derived = ["pop_2", "pop_3"],
    ancestral = "ancestral_2_3"
)

demography.add_population_split(
    time = split_time_1_2_3,
    derived = ["ancestral_2_3", "pop_1"],
    ancestral = "ancestral_1_2_3"
)

demography.add_population_split(
    time = split_time_0_1_2_3,
    derived = ["ancestral_1_2_3", "pop_0"],
    ancestral = "ancestral_0_1_2_3"
)

demography.add_population_split(
    time = split_time_5_6,
    derived = ["pop_5", "pop_6"],
    ancestral = "ancestral_5_6"
)

demography.add_population_split(
    time = split_time_7_5_6,
    derived = ["ancestral_5_6", "pop_7"],
    ancestral = "ancestral_7_5_6"
)

demography.add_population_split(
    time = split_time_8_7_5_6,
    derived = ["ancestral_7_5_6", "pop_8"],
    ancestral = "ancestral_8_7_5_6"
)

demography.add_population_split(
    time = split_time_ingroups,
    derived = ["ancestral_8_7_5_6", "ancestral_0_1_2_3"],
    ancestral = "ancestral_ingroups"
)

demography.add_population_split(
    time = split_time_ingroups_out_2,
    derived = ["ancestral_ingroups", "outpop_1"],
    ancestral = "ancestral_ingroups_out_1"
)

demography.add_population_split(
    time = split_time_ancestral_all,
    derived = ["ancestral_ingroups_out_1", "outpop_0"],
    ancestral = "ancestral_all"
)

# ** Create the admixed population (pop_4)

demography.add_admixture(
    time = admixture_time,
    derived = "pop_4",
    ancestral = ["pop_3", "pop_5"],
    proportions = [0.5, 0.5]
)

# ** Create migration

migration_time = 300
migration_intensity = argue.intensity


demography.add_mass_migration(migration_time, "ancestral_0_1_2_3" , "ancestral_8_7_5_6", )

# ** Sampling

sampling_scheme_0 = msprime.SampleSet( 5, population = "pop_0", time = 0 )
sampling_scheme_1 = msprime.SampleSet( 5, population = "pop_1", time = 0 )
sampling_scheme_2 = msprime.SampleSet( 5, population = "pop_2", time = 0 )
sampling_scheme_3 = msprime.SampleSet( 5, population = "pop_3", time = 0 )
sampling_scheme_4 = msprime.SampleSet( 5, population = "pop_4", time = 0 )
sampling_scheme_5 = msprime.SampleSet( 5, population = "pop_5", time = 0 )
sampling_scheme_6 = msprime.SampleSet( 5, population = "pop_6", time = 0 )
sampling_scheme_7 = msprime.SampleSet( 5, population = "pop_7", time = 0 )
sampling_scheme_8 = msprime.SampleSet( 5, population = "pop_8", time = 0 )

sampling_scheme_out_1 = msprime.SampleSet( 10, population = "outpop_1", time = 0 )
sampling_scheme_out_0 = msprime.SampleSet( 10, population = "outpop_0", time = 0 )

# Merge all the lists of samples
all_samples = [sampling_scheme_0, sampling_scheme_1, sampling_scheme_1, sampling_scheme_3, sampling_scheme_4, sampling_scheme_5, sampling_scheme_6, sampling_scheme_7, sampling_scheme_8, sampling_scheme_out_1, sampling_scheme_out_0]

# ** Sorting events

## sort.events() is mandatoty because they are not passed in chronological order
demography.sort_events()
# print(demography.debug())

# * Output naming

names = [ pathlib.Path( argue.out_folder+"/"+argue.name+"rep_"+str(i)+".vcf" ) for i in range(argue.how_many) ]
names = [ argue.out_folder+"/"+argue.name+"_rep_"+str(i)+".vcf" for i in range(argue.how_many) ]
##  Change "output.vcf"

# * Start the simulation

for i in range(argue.how_many):
    
    ts = msprime.sim_ancestry(
        demography = demography,
        samples = all_samples,
        recombination_rate = 1e-8,
        sequence_length = 5e+6,
        ploidy = 1
    )
    mutated = msprime.sim_mutations( ts, rate = 1e-8, model = 'binary', discrete_genome = True, keep = False )
# ** Modifing individual names for PLiNK integration
    n_dip_indv = int(ts.num_samples)
    indv_names = [f"tsk_{i}indv" for i in range(n_dip_indv)]
    with open(names[i], "w") as vcf_file:
        ts.write_vcf(vcf_file, individual_names=indv_names)
# Plink doesn't like when individuals names end with "_0".
# The previous lines modify sample names to avoid this problem.

