from os import path
from norm.rosen import normalize, normalize_consecutive, normalize_batch
import scipy.io as sio
from scipy.sparse import coo_matrix
from pandas import Series, DataFrame
import matplotlib.pyplot as plt
import seaborn as sns

PRJ_DIR = "/home/derek/research/Kim-Lab/normalization-simulation/"
OUT_DIR = f"{PRJ_DIR}exp/exp-14/out/"
WRT_FLG = True

file_name = f"{OUT_DIR}s.cd4.cd8.eff.mm.mtx"
counts_sparse = sio.mmread(file_name)

counts = counts_sparse.toarray().astype('float64')
res = normalize_consecutive(counts)

# res = normalize_consecutive(counts, verbose=True)
counts_norm = res[0]
slope_list = res[1]
intercept_list = res[2]

if WRT_FLG:
    file_name = f"{OUT_DIR}s.norm.cd4.cd8.eff.mm"
    if not path.exists(file_name):
        print(f'writing to: {file_name}')
        sio.mmwrite(target=file_name, a=coo_matrix(counts_norm))
        # DataFrame(counts_norm).to_csv(file_name, index=False, header=False)

