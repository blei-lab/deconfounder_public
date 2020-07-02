import numpy as np
import pandas as pd

hapmap_gene = pd.read_csv("raw_csv/hapmap_mimno_genes.csv", header=None)
hapmap_pop = pd.read_csv("raw_csv/hapmap_mimno_pops.csv", header=None)

# take out snp values that are larger than 3
hapmap_gene_clean = hapmap_gene.iloc[:,np.where(np.array(hapmap_gene).max(axis=0)<3)[0]]
# take out snps whose values = 2 all the time
hapmap_gene_clean = hapmap_gene_clean.iloc[:,np.where(np.array(hapmap_gene_clean).mean(axis=0)<1.8)[0]]
# original cutoff is 1.9

hapmap_gene_clean.to_csv("clean_csv/genes.csv", index=False)

n_hapmapgenes = hapmap_gene_clean.shape[1]

ps = np.array(hapmap_gene_clean).mean(axis=0) / 2.
np.savetxt("clean_csv/ps.csv", ps, delimiter=",")

# 1-49 ASW 50-161 CEU 162-211 MEX 212-324 YRI
mean_1 = hapmap_gene_clean.iloc[:49, :].mean(axis=0)
mean_2 = hapmap_gene_clean.iloc[49:161, :].mean(axis=0)
mean_3 = hapmap_gene_clean.iloc[161:211, :].mean(axis=0)
mean_4 = hapmap_gene_clean.iloc[211:324, :].mean(axis=0)
mean_all = hapmap_gene_clean.mean(axis=0)

btwngp_var = (49 * (mean_1 - mean_all)**2 + 112 * (mean_2 - mean_all)**2 + 50 * (mean_3 - mean_all)**2 + 113 * (mean_4 - mean_all)**2)
tt_var = hapmap_gene_clean.var(axis=0) * 324

Fs = btwngp_var / tt_var

Fs = np.array(Fs)
np.savetxt("clean_csv/Fs.csv", Fs, delimiter=",")
