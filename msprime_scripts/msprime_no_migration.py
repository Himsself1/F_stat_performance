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

split_time_2_3
split_time_1_2_3
split_time_0_1_2_3
split_time_5_6
split_time_7_5_6
split_time_8_7_5_6
split_time_ingroups
split_time_ingroups_out_2
split_time_ancestral_all


demography = msprime.Demography()

# ** Create populations 0-9

[demography.add_population(name= "pop_"+str(i) ) for i in range(9)]

# ** Create 2 outgroups for reference

[demography.add_population(name= "outpop_"+str(i) ) for i in range(2)]

# ** Create ancestral populations

demography.add_population(name = "ancestral_2_3")
demography.add_population(name = "ancestral_1_2_3")
demography.add_population(name = "ancestral_0_1_2_3")
demography.add_population(name = "ancestral_5_6")
demography.add_population(name = "ancestral_7_5_6")
demography.add_population(name = "ancestral_8_7_5_6")
demography.add_population(name = "ancestral_ingroups")
demography.add_population(name = "ancestral_all")


# ** Create coalescent events

demography.add_population_split(
    time = split_time_2_3,
    derived = [ "pop_2", "pop_3"],
    ancestral = "ancestral_2_3"
)

demography.add_population_split(
    time = split_time_1_2_3,
    derived = [ "ancestral_2_3", "pop_1"],
    ancestral = "ancestral_1_2_3"
)

demography.add_population_split(
    time = split_time_0_1_2_3,
    derived = [ "ancestral_1_2_3", "pop_0"],
    ancestral = "ancestral_0_1_2_3"
)

demography.add_population_split(
    time = split_time_5_6,
    derived = [ "pop_5", "pop_6"],
    ancestral = "ancestral_5_6"
)

demography.add_population_split(
    time = split_time_7_5_6,
    derived = [ "ancestral_5_6", "pop_7"],
    ancestral = "ancestral_7_5_6"
)

demography.add_population_split(
    time = split_time_8_7_5_6,
    derived = [ "ancestral_7_5_6", "pop_8"],
    ancestral = "ancestral_8_7_5_6"
)

demography.add_population_split(
    time = split_time_ingroups,
    derived = [ "ancestral_8_7_5_6", "ancestral_0_1_2_3"],
    ancestral = "ancestral_ingroups"
)

demography.add_population_split(
    time = split_time_ingroups_out_2,
    derived = [ "ancestral_ingroups", "outpop_2"],
    ancestral = "ancestral_outpop_2"
)

demography.add_population_split(
    time = split_time_ancestral_all,
    derived = [ "ancestral_outpop_2", "outpop_1"],
    ancestral = "ancestral_all"
)

# ** Create the admixed population (pop_4)
