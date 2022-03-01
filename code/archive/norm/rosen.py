import numpy as np
from scipy import stats
import pandas as pd
from pandas import Series, DataFrame

from norm.resit_fit import linear_resist_fit_p, linear_resist_fit_robust, linear_resist_fit_robust_mix


def get_transcript_centroid(counts):
    return np.mean(counts, axis=1)


# We have to add regularization entry or the slope will eventually be 0
def normalize(counts, epsilon=0.001, verbose=False, p_default=0.5):
    n_gene, n_cell = counts.shape
    print(f'n_gene: {n_gene}')
    print(f'n_cell: {n_cell}')

    t_c = counts[:, 0]
    # delta = |d' - d|
    d = 0
    d_prime = np.inf
    delta = np.inf

    counter = 0

    norm_counts = counts.copy()

    while delta > epsilon:
        i = counter % n_cell
        # norm_counts[:, i] = linear_resist_fit_p(t_c, norm_counts[:, i], p)
        # norm_counts[:, i] = linear_resist_fit_p(norm_counts[:, i], t_c, p)
        # norm_counts[:, i] = linear_resist_fit_robust(norm_counts[:, i], t_c, p)
        norm_counts[:, i] = linear_resist_fit_robust_mix(norm_counts[:, i], t_c, p_default=p_default, verbose=verbose)
        # norm_counts[:, i] = linear_resist_fit_robust_mix(norm_counts[:, i], t_c, verbose=True)

        t_c_prime = get_transcript_centroid(norm_counts)
        d = np.linalg.norm(t_c_prime - t_c)
        delta = abs(d_prime - d)
        # Update parameters
        d_prime = d
        t_c = t_c_prime
        counter = counter + 1

        # if counter % 100 == 0:
        print(f"counter: {counter}; d: {d}, delta: {delta}")
    return norm_counts


def normalize_consecutive(counts, p_default=0.5, verbose=False):
    n_gene, n_cell = counts.shape
    n_iter = 1
    print(f'n_gene: {n_gene}')
    print(f'n_cell: {n_cell}')

    counter = 0
    norm_counts = counts.copy()
    slope_list = []
    intercept_list = []
    for j in range(n_iter):
        print(f'j: {j}')
        t_c = get_transcript_centroid(norm_counts)
        print(f'seq-depth of tc: {np.sum(t_c)}')
        for i in range(n_cell):
            if verbose:
                print(f'---------{i}---------')
            else:
                if counter % 100 == 0:
                    print(f'i: {i}')
            res = linear_resist_fit_robust_mix(norm_counts[:, i], t_c, p_default=p_default, verbose=verbose)

            norm_counts[:, i] = res[0]
            slope_list.append(res[1])
            intercept_list.append(res[2])

            counter = counter + 1
            # print(f'# of 0s: {np.sum(counts[:, i]==0)}: {np.sum(norm_counts[:, i]==0)}')
    return norm_counts, slope_list, intercept_list


def normalize_batch(counts, p_default=0.5, verbose=False, batch=None):
    n_gene, n_cell = counts.shape
    n_iter = 3
    print(f'n_gene: {n_gene}')
    print(f'n_cell: {n_cell}')
    print(f'n_iter: {n_iter}')

    batch = np.array(batch)
    first_type = np.unique(batch)[0]
    first_type_idx = (batch == first_type)
    
    counter = 0
    norm_counts = counts.copy()
    slope_list = []
    intercept_list = []
    for j in range(n_iter):
        print(f'j: {j}')
        t_c = get_transcript_centroid(norm_counts[:, first_type_idx])
        print(f'seq-depth of tc: {np.sum(t_c)}')
        for i in range(n_cell):
            if verbose:
                print(f'---------{i}---------')
            else:
                if counter % 100 == 0:
                    print(f'i: {i}')
            res = linear_resist_fit_robust_mix(norm_counts[:, i], t_c, p_default=p_default, verbose=verbose)

            norm_counts[:, i] = res[0]
            slope_list.append(res[1])
            intercept_list.append(res[2])

            counter = counter + 1
            # print(f'# of 0s: {np.sum(counts[:, i]==0)}: {np.sum(norm_counts[:, i]==0)}')
    return norm_counts, slope_list, intercept_list

def main():
    file_name = "/home/derek/research/Kim-Lab/normalization-simulation/exp/exp-7/out/sim-2/true_counts.csv"
    counts = pd.read_csv(file_name, index_col=0).to_numpy()
    # print(f"Estimated p:{p}")
    x = normalize(counts)

    file_name = "/home/derek/research/Kim-Lab/normalization-simulation/exp/exp-7/out/sim-2/norm_counts.csv"
    DataFrame(x).to_csv(file_name, index=False, header=False)


if __name__ == '__main__':
    main()
