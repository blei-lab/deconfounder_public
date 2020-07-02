import numpy as np
import fastStructure 
import parse_bed
import parse_str
import random
import sys
import warnings
import pandas as pd

from numpy.random import dirichlet, beta

import itertools

params = {'inputfile': 'mimno_hapmap3/hapmap3',
        }

G = parse_bed.load(params['inputfile'])
G = np.require(G, dtype=np.uint8, requirements='C')

# G = G[:,0:1000]
G.shape

np.savetxt("hapmap_mimno_genes.csv", G, delimiter=",")


pops = pd.read_csv('mimno_hapmap3/hapmap3.pops', header=None)
pops = pops.values
np.savetxt("hapmap_mimno_pops.csv", pops, delimiter=",")