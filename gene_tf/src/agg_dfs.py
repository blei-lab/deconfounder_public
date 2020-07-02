import os
from fnmatch import fnmatch
import pandas as pd
import numpy as np

root='./gene_tf/res/causalest'

for simset in ["BN", "TGP", "HGDP", "PSD", "SP"]:
    for alpha in [1, 10, 25, 50, 100]:
        for n_units in [5000]:
            for n_causes in [5000, 100000]:
                for snpsig in [10, 40]:
                    for confint in [20, 40]:
                        for causalprop in [10]:
                            filenames = []
                            if simset in ["BN", "TGP", "HGDP"]:
                                pattern = simset + '_units'+str(n_units) + '_causes'+str(n_causes) + "*" + "_cv0" + "_snpsig" + str(snpsig) + "_confint" + str(confint) + "_cp" + str(causalprop) + "*.csv"
                            elif simset in ["PSD", "SP"]:
                                pattern = simset + 'alpha'+str(alpha) + '_units'+str(n_units) + '_causes'+str(n_causes) + "*" + "_cv0" + "_snpsig" + str(snpsig) + "_confint" + str(confint) + "_cp" + str(causalprop) + "*.csv"
                            print(pattern)

                            for path, subdirs, files in os.walk(root):
                                for name in files:
                                    if fnmatch(name, pattern):
                                        filenames.append(os.path.join(path, name)) 
                            print(filenames)

                            dataframes = [pd.read_csv(df) for df in filenames]

                            print(len(dataframes), filenames)

                            if len(dataframes) > 0:
                                df_concat = pd.concat(dataframes)
                                result = df_concat.groupby('adjustments').mean().sort_values('Unnamed: 0').iloc[:,1:]
                                if simset in ["BN", "TGP", "HGDP"]:
                                    result_fn = simset + '_units'+str(n_units) + '_causes'+str(n_causes) + "_cv0" + "_snipsig" + str(snpsig) + "_confint" + str(confint) + "_cp" + str(causalprop) + "_avg.csv"
                                elif simset in ["PSD", "SP"]:
                                    result_fn = simset + 'alpha'+str(alpha) + '_units'+str(n_units) + '_causes'+str(n_causes) + "_cv0" + "_snipsig" + str(snpsig) + "_confint" + str(confint) + "_cp" + str(causalprop) + "_avg.csv"
                                result.to_csv(root + "/" + result_fn)

                                

