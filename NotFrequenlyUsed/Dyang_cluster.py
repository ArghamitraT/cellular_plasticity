import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib as mpl
params = {
    'font.family': "Helvetica",
    'figure.dpi': 300
   }
mpl.rcParams.update(params)
mpl.rc('savefig', dpi=300)

data_directory = "/path/to/KPTracer-Data"

adata = sc.read_h5ad(f"{data_directory}/expression/adata_processed.nt.h5ad")

sigscores = pd.read_csv(f"../Figure6_S6/data/fitness_signature_scores.tsv", sep='\t', index_col = 0, usecols=['FitnessSignature_NT'])