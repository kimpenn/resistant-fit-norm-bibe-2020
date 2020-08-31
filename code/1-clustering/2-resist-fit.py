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

file_name = f"{OUT_DIR}cm.csv"
counts_sparse = sio.mmread(file_name)
counts = counts_sparse.toarray().astype('float64')

# counts = counts_sparse.todense()
# counts = pd.read_csv(file_name, index_col=0)
# counts = counts.to_numpy(dtype=np.float64)

res = normalize_consecutive(counts)
# res = normalize_consecutive(counts, verbose=True)
counts_norm = res[0]
slope_list = res[1]
intercept_list = res[2]

if WRT_FLG:
    # file_name = f"{OUT_DIR}norm.cm-iter-1-no-intercept.mm"
    file_name = f"{OUT_DIR}norm.cm-iter-1-intercept.mm"
    if not path.exists(file_name):
        print(f'writing to: {file_name}')
        sio.mmwrite(target=file_name, a=coo_matrix(counts_norm))
        # DataFrame(counts_norm).to_csv(file_name, index=False, header=False)

group_list = ['group1']*2776 + ['group2']*2775
slope_df = DataFrame({
  'slope': slope_list,
  'group': group_list
})
ax = sns.boxplot(x='group', y='slope', data=slope_df, orient='v')
plt.show()
# counts_norm = normalize(counts, verbose=True)
ax = sns.boxplot(x=slope_list[0:2776], orient='v')
plt.show()
ax = sns.boxplot(x=slope_list[2776:5551], orient='v')
plt.show()
ax = sns.boxplot(x=intercept_list, orient='v')
plt.show()


# def get_transcript_centroid(counts):
#     return np.mean(counts, axis=1)
#
# # -----------------------------
# # Pick one cell from monocytes
# # -----------------------------
# t_c = get_transcript_centroid(counts)
# cell_select_1 = counts[:, 0].copy()
# cell_norm_1 = linear_resist_fit_robust_mix(cell_select_1, t_c, verbose=True)
#
# plt.subplot()
# plt.plot(cell_select_1, t_c, 'o')
# plt.show()
# plt.subplot()
# plt.plot(cell_norm_1, t_c, 'o')
# plt.show()
# # -----------------------------
# # Pick one cell from down-sampled monocytes
# # -----------------------------
# cell_select_2 = counts[:, 3000].copy()
# cell_norm_2 = linear_resist_fit_robust_mix(cell_select_2, t_c, verbose=True)
#
# plt.subplot()
# plt.plot(cell_select_2, t_c, 'o')
# plt.show()
#
# plt.subplot()
# plt.plot(cell_norm_2, t_c, 'o')
# plt.show()
# # -----------------------------
# # Bioplot for two selected cells
# # -----------------------------
# plt.subplot()
# plt.plot(cell_norm_1, cell_norm_2, 'o')
# plt.show()
#
# # -----------------------------
# # Bioplot for two selected cells
# # -----------------------------
# plt.subplot()
# plt.plot(counts_norm[:, 0], counts_norm[:, 3000], 'o')
# plt.show()
#
# plt.subplot()
# plt.plot(counts_norm[:, 42], counts_norm[:, 3099], 'o')
# plt.show()
#
# plt.subplot()
# plt.plot(counts_norm[:, 42], counts_norm[:, 3588], 'o')
# plt.show()