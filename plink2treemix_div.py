import sys
import pandas as pd

if len(sys.argv) != 2:
    sys.stderr.write("\n")
    sys.stderr.write(" (script will write uncompressed output to stdout)\n")
    sys.stderr.write("\n")
    sys.exit()

plink_freq_file = sys.argv[1]

freq_file = "/home/thom_nelson/LewCard/seq/popgenetics/mapped2CE10chromonome/plink/diversity_GATK_CE10chromonome.8.variable.SNPs.frq.strat_toy"

freq_df = pd.read_csv(plink_freq_file, delim_whitespace=True)

snp_IDs = freq_df["SNP"].unique()
pop_IDs = freq_df["CLST"].unique()
nsnps = len(snp_IDs)
npops = len(pop_IDs)

sys.stderr.write("SNPs in file: %s\n" % nsnps)
sys.stderr.write("Populations in input (n = %s):\n" % npops)
[sys.stderr.write("  %s\n" % pop) for pop in pop_IDs]

pop_IDs_out = " ".join(pop_IDs)
sys.stdout.write("%s\n" % pop_IDs_out)

# record the snp order
snp_order = freq_df.drop_duplicates("SNP")["SNP"].tolist()

# transform MAC & NCHROBS by pivot
mac_wide = freq_df.pivot(index="SNP", columns="CLST", values="MAC")
nobs_wide = freq_df.pivot(index="SNP", columns="CLST", values="NCHROBS")

# reorder
mac_wide = mac_wide[pop_IDs]
nobs_wide = nobs_wide[pop_IDs]

# mask NaN
mask = mac_wide.notna().all(axis=1) & nobs_wide.notna().all(axis=1)
mac_wide = mac_wide[mask]
nobs_wide = nobs_wide[mask]

# reorder
snp_order2 = [s for s in snp_order if s in mac_wide.index]
mac_wide  = mac_wide.loc[snp_order2]
nobs_wide = nobs_wide.loc[snp_order2]

# major
minor = mac_wide.to_numpy().astype("int32")
nalleles = nobs_wide.to_numpy().astype("int32")
major = nalleles - minor

# output
nsnps_out = 0
for maj_row, min_row in zip(major, minor):
    allele_counts_per_pop = [
        "{},{}".format(maj, min_) for maj, min_ in zip(maj_row, min_row)
    ]
    sys.stdout.write(" ".join(allele_counts_per_pop) + "\n")
    nsnps_out += 1
    if nsnps_out % 100000 == 0:
        sys.stderr.write("SNPs written: {}\r".format(nsnps_out))

sys.stderr.write("\nTotal SNPs written (chr12): {}\n".format(nsnps_out))
