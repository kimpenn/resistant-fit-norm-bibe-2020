import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import matplotlib.pyplot as plt
from os import path
import seaborn as sns
from norm.resit_fit import linear_resist_fit_robust_mix, linear_resist_fit_p
from norm.rosen import normalize, normalize_consecutive
import scipy.io as sio
from scipy.sparse import coo_matrix

PRJ_DIR = "/home/derek/research/Kim-Lab/normalization-simulation/"
OUT_DIR = f"{PRJ_DIR}exp/exp-12/out/"
WRT_FLG = True

file_list = ["CD4_Memory", "CD4_Naive", "Mono_CD14", "B_Pre", "CD8_Memory",
             "CD8_Naive", "NK_Dim", "Mono_CD16", "CD8_Effector", "B_Pro",
             "DC", "NK_Bright", "Mk", "pDC"]

for file_prefix in file_list:
    print(f'working on {file_prefix}......')
    # Read downsampled count matrix
    file_name = f"{OUT_DIR}{file_prefix}-cm.mm"
    counts_sparse = sio.mmread(file_name)
    counts = counts_sparse.toarray().astype('float64')
    res = normalize_consecutive(counts)
    
    counts_norm = res[0]
    slope_list = res[1]
    intercept_list = res[2]

    if WRT_FLG:
        file_name = f"{OUT_DIR}{file_prefix}cm-norm.mm"
        if not path.exists(file_name):
            print(f'writing to: {file_name}')
            sio.mmwrite(target=file_name, a=coo_matrix(counts_norm))
            # DataFrame(counts_norm).to_csv(file_name, index=False, header=False)
