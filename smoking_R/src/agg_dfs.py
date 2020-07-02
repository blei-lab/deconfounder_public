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

for factor_model in ["quadratic", "linear"]:
    for factor_model_K in [1, 2, 3]:
        for outcome_model in ["linear"]:
            for simset in ["indep", "dep"]:
                for estctrl in ["Zctrl", "Actrl", "Acovctrl", "Zctrl", "nctrl", "oracle"]:
                    for dataseed in ["20133224"]:
                        for factorseed in ["20144226", "20144133", "20144042", "20144219", "20144148", "20135747", "20144211", "20144201", "20144119"]:
                            for randseed in ['20193837', '20235400', '21061444', '21061558', '21111450',
                               '21155118', '21183010', '21195346', '22094716', '22121550',
                               '22121631', '22121649', '22121713', '22121730', '22145351',
                               '22145437', '22145457', '22145515', '22145533', '22155052',
                               '22180017', '22195956', '22212357', '23002709']: # this is estimation seed
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

                                dataframes = [pd.read_csv(df) for df in filenames]

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



