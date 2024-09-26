# * Libraries

import msprime
from pathlib import Path
import argparse

# * Command line arguments

cli = argparse.ArgumentParser(
    description = 'msprime simulations with ',
    prog = 'simulations_msprime.py'
)

# * Build population model 

# ** Configure split times

## Need to configure these split times in order to finalize the model
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

pop_size_0 = 1000
pop_size_1 = 1000
pop_size_2 = 1000
pop_size_3 = 1000
pop_size_4 = 1000
pop_size_5 = 1000
pop_size_6 = 1000
pop_size_7 = 1000
pop_size_8 = 1000
pop_size_out_1 = 1000
pop_size_out_0 = 1000

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

[demography.add_population(name= "pop_"+str(i) ) for i in range(9)]

# ** Create 2 outgroups for reference

[demography.add_population(name= "outpop_"+str(i) ) for i in range(2)]

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
    ancestral = "ancestral_outpop_2"
)

demography.add_population_split(
    time = split_time_ancestral_all,
    derived = ["ancestral_outpop_1", "outpop_0"],
    ancestral = "ancestral_all"
)

# ** Create the admixed population (pop_4)

demography.add_admixture(
    time = admixture_time,
    derived = "pop_4",
    ancestral = ["pop_3", "pop_5"],
    proportions = [0.5, 0.5]
