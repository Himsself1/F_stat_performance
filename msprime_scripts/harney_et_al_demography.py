# * Libraries

import msprime
import argparse
import pathlib
import numpy
import timeit
import os
import sys
import importlib.util

# * Command line arguments

cli = argparse.ArgumentParser(
    description = 'msprime simulations with model my Harney et al. 2021',
    prog = 'harney_et_al_demography.py'
)

cli.add_argument("-out_folder", type = str, required = True,
                 help = "Full path to output folder. The program will create all child folders that may not exist.")
cli.add_argument("-name", type = str, required = True,
                 help = "Name prefix of output vcf files.")
cli.add_argument("-how_many", type = int, default = 1,
                 help = 'How many simulations are going to be run in one go.')
cli.add_argument("-rec", type = float, default = 1.25e-8,
                 help = 'Recombination rate.')
cli.add_argument("-mut", type = float, default = 1.25e-8,
                 help = 'Mutation rate.')
cli.add_argument("-seq_length", type = float, default = 5e+6,
                 help = 'Genome Size.')
cli.add_argument("-nsamples", type = int, default = 5,
                 help = 'Number of samples per population.')
cli.add_argument("-mix_prop_14", type = float, default = 0.5,
                 help = "Admixture proportion of 14a. Proportion of 14b is 1-`mix_prop_14`")
cli.add_argument("-mix_prop_15", type = float, default = 0.55,
                 help = "Admixture proportion of 15a. Proportion of 14b is 1-`mix_prop_15`")
cli.add_argument("-path_to_custom_module", type = str,
                 default = "/home/stefanos/F_stat_performance/msprime_scripts/ts_to_eigenstrat.py",
                 help = "Custom module for EIGENSTRAT output.")

argue = cli.parse_args()

# * Debug Arguments

# argue = cli.parse_args([
#     "-out_folder", "/home/stefanos/test_python_script",
#     "-name", "test_harney"   
# ])

# ** Load custom module

spec = importlib.util.spec_from_file_location("ts_to_eigenstrat", pathlib.Path(argue.path_to_custom_module))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
ts_to_eigenstrat = module.ts_to_eigenstrat

# ** Build Directory Tree

run_folder = pathlib.Path(argue.out_folder, argue.name)
eigenstrat_folder = pathlib.Path(run_folder, "eig")
pathlib.Path(run_folder).mkdir(parents = True, exist_ok = True) # Makes Output directories
pathlib.Path(eigenstrat_folder).mkdir(parents = True, exist_ok = True) # Makes Output directories

# * Build population model 

# ** Population Sizes
## Kept the same asignment scheme as Harney for debugging purposes.

## The only difference is that the legacy msprime.simulate() which is
## originally used simulates monoploid samples. Thus, I need to halve
## population sizes
N_A        = 50000
N_B        = 12000
N_Bp       = 10000
N_BBp      = 12000
N_C        = 25000
N_Cp       = 10000
N_X        = 10000
N_CCp      = 3300
N_O        = 80000
N_ABBp     = 7000
N_ABBpCCp  = 2500
N_ABBpCCpO = 10000

N_0 = N_A 
N_1 = N_A
N_2 = N_A
N_3 = N_A
N_4 = N_A
N_5 = N_B 
N_6 = N_B
N_7 = N_ABBp
N_8 = N_CCp
N_9 = N_C 
N_10 = N_ABBpCCp
N_11 = N_ABBpCCp
N_12 = N_ABBpCCp
N_13 = N_O 
N_14 = N_X 
N_15 = N_ABBpCCp
N_16 = 1 
N_17 = 1

list_of_pop_sizes = [N_0, N_1, N_2, N_3, N_4, N_5, N_6, N_7, N_8, N_9, N_10, N_11, N_12, N_13, N_14, N_15]

demography = msprime.Demography()

# ** Populations 0-16

[demography.add_population(name = "pop_"+str(i), initial_size = list_of_pop_sizes[i]) for i in range(16)]
# Populations 14a/b and 15a/b will not get any present-day samples.
demography.add_population(name = "pop_14a", initial_size = N_14, initially_active = False)
demography.add_population(name = "pop_14b", initial_size = N_14, initially_active = False)
demography.add_population(name = "pop_15a", initial_size = N_15, initially_active = False)
demography.add_population(name = "pop_15b", initial_size = N_15, initially_active = False)

# ** Ancsetral Populations

# *** Leaf Populations

demography.add_population(name = "anc_0_1", initial_size = N_A)
demography.add_population(name = "anc_0_1_2", initial_size = N_A)
demography.add_population(name = "anc_0_1_2_3", initial_size = N_A)
demography.add_population(name = "anc_0_1_2_3_4", initial_size = N_A)

demography.add_population(name = "anc_5_6", initial_size = N_B)  

demography.add_population(name = "anc_9_14b", initial_size = N_CCp)

demography.add_population(name = "anc_11_12", initial_size = N_ABBpCCp)
demography.add_population(name = "anc_11_12_15b", initial_size = N_ABBpCCp) # ((11,12),15b)

demography.add_population(name = "anc_0_1_2_3_4_anc_5_6_14a_7_anc_8_9_14b_15a",
                          initial_size = N_ABBpCCp) # ((((((((0,1),2),3),4),((5,6),14a)),7),((9,14a),8)),15a);

demography.add_population(name = "anc_0_1_2_3_4_anc_5_6_14a_7_anc_8_9_14b_15a_10",
                          initial_size = N_ABBpCCp) # (((((((((0,1),2),3),4),((5,6),14a)),7),((9,14a),8)),15a),10);

demography.add_population(name = "anc_all", initial_size = N_ABBpCCpO) # (((((((((((0,1),2),3),4),((5,6),14a)),7),((9,14a),8)),15a),10),((11,12),15b)),13);

# *** Internal Populations

demography.add_population(name = "anc_5_6_14a", initial_size = N_BBp)
demography.add_population(name = "anc_0_1_2_3_4_anc_5_6_14a", initial_size = N_ABBp) # (((((0,1),2),3),4),(5,6),14a);
demography.add_population(name = "anc_0_1_2_3_4_anc_5_6_14a_7", initial_size = N_ABBp) # ((((((0,1),2),3),4),((5,6),14a)),7);

demography.add_population(name = "anc_8_9_14b", initial_size = N_CCp)

demography.add_population(name = "anc_0_1_2_3_4_anc_5_6_14a_7_anc_8_9_14b",
                          initial_size = N_ABBpCCp) # (((((((0,1),2),3),4),((5,6),14a)),7),((9,14a),8));

demography.add_population(name = "anc_0_1_2_3_4_anc_5_6_14a_7_anc_8_9_14b_15a_10_anc_11_12_15b",
                          initial_size = N_ABBpCCp) # ((((((((((0,1),2),3),4),((5,6),14a)),7),((9,14b),8)),15a),10),((11,12),15b));

# ** Configure Splits and Split Times

# *** Population Splits

T_0_13 = 2400
T_0_11 = 2000
T_0_10 = 1600
T_0_8 = 1200
T_11_12	= 1000
T_0_7 = 800
T_8_9 = 600
T_0_5 = 400
T_0_4 = 300
T_0_3 = 200
T_5_6 = 150
T_0_2 =  100
T_0_1 = 50

# *** Split times of Sources of 14 Admixture
T_5_14a = 240
T_9_14b = 280


# *** Split times of Sources of 15 Admixture
T_0_15a = 1450
T_11_15b = 1450

demography.add_population_split(time = T_0_1, derived = ["pop_0", "pop_1"], ancestral = "anc_0_1") # (0,1);
demography.add_population_split(time = T_0_2, derived = ["anc_0_1", "pop_2"], ancestral = "anc_0_1_2") # ((0,1),2);
demography.add_population_split(time = T_0_3, derived = ["anc_0_1_2", "pop_3"], ancestral = "anc_0_1_2_3") # (((0,1),2),3);
demography.add_population_split(time = T_0_4, derived = ["anc_0_1_2_3", "pop_4"],
                                ancestral = "anc_0_1_2_3_4")  # ((((0,1),2),3),4);

demography.add_population_split(time = T_5_6, derived = ["pop_5", "pop_6"], ancestral = "anc_5_6") # (5,6);
demography.add_population_split(time = T_5_14a, derived = ["anc_5_6", "pop_14a"], ancestral = "anc_5_6_14a") # ((5,6),14a);
demography.add_population_split(time = T_0_5, derived = ["anc_0_1_2_3_4", "anc_5_6_14a"],
                                ancestral = "anc_0_1_2_3_4_anc_5_6_14a") # (((((0,1),2),3),4),(5,6),14a);
demography.add_population_split(time = T_0_7, derived = ["anc_0_1_2_3_4_anc_5_6_14a", "pop_7"],
                                ancestral = "anc_0_1_2_3_4_anc_5_6_14a_7") # ((((((0,1),2),3),4),((5,6),14a)),7);

demography.add_population_split(time = T_9_14b, derived = ["pop_9", "pop_14b"], ancestral = "anc_9_14b") # (9,14b);
demography.add_population_split(time = T_8_9, derived = ["anc_9_14b", "pop_8"], ancestral = "anc_8_9_14b") # ((9,14b),8);

demography.add_population_split(time = T_0_8,
                                derived = ["anc_0_1_2_3_4_anc_5_6_14a_7", "anc_8_9_14b"],
                                ancestral = "anc_0_1_2_3_4_anc_5_6_14a_7_anc_8_9_14b") # (((((((0,1),2),3),4),((5,6),14a)),7),((9,14a),8));

demography.add_population_split(time = T_0_15a,
                                derived = ["anc_0_1_2_3_4_anc_5_6_14a_7_anc_8_9_14b", "pop_15a"],
                                ancestral = "anc_0_1_2_3_4_anc_5_6_14a_7_anc_8_9_14b_15a") # ((((((((0,1),2),3),4),((5,6),14a)),7),((9,14a),8)),15a);

demography.add_population_split(time = T_0_10,
                                derived = ["anc_0_1_2_3_4_anc_5_6_14a_7_anc_8_9_14b_15a", "pop_10"],
                                ancestral = "anc_0_1_2_3_4_anc_5_6_14a_7_anc_8_9_14b_15a_10") # (((((((((0,1),2),3),4),((5,6),14a)),7),((9,14a),8)),15a),10);

demography.add_population_split(time = T_11_12, derived = ["pop_11", "pop_12"], ancestral = "anc_11_12") # (0,1);
demography.add_population_split(time = T_11_15b, derived = ["anc_11_12", "pop_15b"],
                                ancestral = "anc_11_12_15b") # ((11,12),15b);

demography.add_population_split(time = T_0_11,
                                derived = ["anc_0_1_2_3_4_anc_5_6_14a_7_anc_8_9_14b_15a_10", "anc_11_12_15b"],
                                ancestral = "anc_0_1_2_3_4_anc_5_6_14a_7_anc_8_9_14b_15a_10_anc_11_12_15b") # ((((((((((0,1),2),3),4),((5,6),14a)),7),((9,14a),8)),15a),10),((11,12),15b));

demography.add_population_split(time = T_0_13,
                                derived = ["anc_0_1_2_3_4_anc_5_6_14a_7_anc_8_9_14b_15a_10_anc_11_12_15b", "pop_13"],
                                ancestral = "anc_all") # (((((((((((0,1),2),3),4),((5,6),14a)),7),((9,14a),8)),15a),10),((11,12),15b)),13);

# ** Admixture Events

# Make proportions a command-line variable.
T_admix_14 = 40
T_admix_15 = 1400

prop_14 = argue.mix_prop_14
prop_15 = argue.mix_prop_15

demography.add_admixture(time = T_admix_14,
                         derived = "pop_14",
                         ancestral = ["pop_14a","pop_14b"],
                         proportions = [argue.mix_prop_14, 1-argue.mix_prop_14])

demography.add_admixture(time = T_admix_15,
                         derived = "pop_15",
                         ancestral = ["pop_15a","pop_15b"],
                         proportions = [argue.mix_prop_15, 1-argue.mix_prop_15])

demography.sort_events()
# print(demography.debug())

# ** Sample Schemes

n_samples = argue.nsamples
all_samples = [msprime.SampleSet(n_samples, population = "pop_{}".format(i)) for i in range(16)]

# * Starting the Simuation

# ** Naming

prefix = pathlib.Path(eigenstrat_folder, argue.name)
ind_metadata_filename = str(pathlib.Path(eigenstrat_folder).parent)+"/"+argue.name+"_model_metadata.tsv"

# ** Loop the simulation

for rep in range(argue.how_many):
    start = timeit.default_timer()
    ts = msprime.sim_ancestry(
        demography = demography,
        samples = all_samples,
        model=[
            # msprime.DiscreteTimeWrightFisher(duration=40),
            msprime.StandardCoalescent(),
        ],
        recombination_rate = argue.rec,
        sequence_length = argue.seq_length,
        ploidy = 2
    )
    stop1 = timeit.default_timer()
    # print('coalescent time: ', stop1 - start)
    mutated = msprime.sim_mutations(ts, rate = argue.mut, model = 'binary', discrete_genome = True, keep = False)
    stop2 = timeit.default_timer()
    ts_to_eigenstrat(
        ts = mutated,
        demography = demography,
        out_prefix = f"{prefix}_rep_{rep}"
    )

# ** Output Sample Metadata for the simulation

n_dip_indv = int(mutated.num_individuals)
indv_names = [f"tsk_{ind.id}_indv" for ind in mutated.individuals()]

with open(ind_metadata_filename, 'w') as metadata:
    metadata.write("Ind_ID\tPopulation\n")
    for i, ind_id in enumerate(mutated.individuals()):
        metadata.write(f"{indv_names[i]}\t{demography.populations[ind_id.population].name}\n")
metadata.close()
