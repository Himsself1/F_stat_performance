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

# ** Create coalescent events

demography.add_population_split(
    time = split_time_2_3,
    derived = [ "pop_2", "pop_3"],
    ancestral=ancestral = "ancestral_2_3"
)

# ** Create the admixed population (pop_4)
