import tskit
import numpy
import pathlib

def ts_to_eigenstrat(ts, demography, out_prefix):
    """
    Export tree sequence to EIGENSTRAT format (.geno, .snp, .ind)
    """

    samples =  ts.individuals()
    n_samples = len(samples)
    n_snps = ts.num_sites

    # .ind file: sample_id sex population
    with open(f"{out_prefix}.ind", "w") as f:
        for i, sample in enumerate(samples):
            pop_name = demography.populations[sample.population].name
            f.write(f"tsk_{sample.id}_indv\tU\t{pop_name}\n")
    f.close()
    print(f"{out_prefix}.ind completed.")
    # .snp file: snp_id chrom genetic_pos phys_pos ref alt
    with open(f"{out_prefix}.snp", "w") as f:
        for site in ts.sites():
            ancestral = site.ancestral_state
            derived = site.mutations[0].derived_state
            f.write(
                f"snp_{int(site.id)}\t1\t{site.position / 1e8:.6f}\t"
                f"{int(site.position)}\t{ancestral}\t{derived}\n"
            )
    f.close()
    print(f"{out_prefix}.snp completed.")
    # .geno file: one row per SNP, one column per individual (0/1/2/9)
    G = ts.genotype_matrix()  # shape: (n_sites, n_samples)
    G_diploid = G[:,range(0,160,2)] + G[:,range(1,160,2)]
    numpy.savetxt(pathlib.Path(f"{out_prefix}.geno"), G_diploid, fmt="%d", delimiter = "")
    print(f"{out_prefix}.geno completed.")

