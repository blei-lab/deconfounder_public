import os
from fnmatch import fnmatch
import pandas as pd
import numpy as np

# accumulate fit random seed in the /src director
pattern = "*.out"
outfilenames = []
for path, subdirs, files in os.walk('./'):
    for name in files:
        if fnmatch(name, pattern):
            outfilenames.append(os.path.join(path, name)) 

randseeds = np.array([outfilename.split('_')[-5] for outfilename in outfilenames])
uniq_randseeds = np.unique(randseeds)

dataseeds = np.array([int(outfilename.split('_')[1]) for outfilename in outfilenames if outfilename.split('_')[1].isdigit()])
uniq_dataseeds = np.unique(dataseeds)


factorseeds = np.array([int(outfilename.split('_')[-6]) for outfilename in outfilenames if outfilename.split('_')[-6].isdigit()])
uniq_factorseeds = np.unique(factorseeds)

for factor_model in ["quadratic", "linear"]:
    for factor_model_K in [1, 2, 3]:
        for outcome_model in ["linear"]:
            for simset in ["indep", "dep"]:
                for estctrl in ["Zctrl", "Actrl", "Zcovctrl", "Acovctrl", "nctrl", "oracle"]:
                    for dataseed in uniq_dataseeds:
                        for factorseed in uniq_factorseeds:
                            for randseed in uniq_randseeds: # this is estimation seed
                                root = '../res/causalest/' + simset + "/" + \
                                    factor_model + "_" + str(factor_model_K)
                                model = estctrl + "_*_" + str(randseed) \
                                    + "_" + str(dataseed)+"_"+str(factorseed)

                                pattern = model + "_res.csv"
                                print(model)

                                filenames = []

                                for path, subdirs, files in os.walk(root):
                                    for name in files:
                                        if fnmatch(name, pattern):
                                            filenames.append(os.path.join(path, name)) 

                                dataframes = [pd.read_csv(df) for df in filenames if np.array(pd.read_csv(df)<1e10).all()]

                                print(len(dataframes), filenames)

                                if len(dataframes) > 0:
                                    df_concat = pd.concat(dataframes)
                                    result = df_concat.groupby(level=0).mean()
                                    result_fn = estctrl \
                                        + "_" + str(randseed) + "_" + \
                                        str(dataseed)+"_"+str(factorseed) \
                                        + "_avg.csv"
                                    result.to_csv(root + "/" + result_fn)
                                    result.to_csv('../res/causalest/' + \
                                        simset + "_" + factor_model + "_" + \
                                        str(factor_model_K) + "_" + result_fn)

                                    print(model, result)



