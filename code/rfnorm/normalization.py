import logging.config
from typing import Tuple

import numpy as np
from scipy.sparse import csr_matrix
from tqdm import tqdm
from rfnorm._resit_fit import _resistant_fit_linear

logging.config.fileConfig(fname='logger.conf', disable_existing_loggers=False)
logger = logging.getLogger(__name__)


def _get_centroid(counts: np.ndarray, by: str = "median") -> np.ndarray:
    """
    Get the center of the dataset as the reference. All the cells will be aligned to it.
    :param counts: Gene x Cell numpy Matrix
    :param by: the method of getting the centroid
    :return: reference cell: Gene x 1 numpy Matrix
    """
    if by == "median":
        return np.median(counts, axis=1)
    elif by == "mean":
        return np.mean(counts, axis=1)
    elif by == "sum":
        return np.sum(counts, axis=1)
    else:
        raise NameError("by must be chosen from list ['median', 'mean', 'sum']")


def rf_normalize_data(counts: np.ndarray,
                      p: np.float32 = 0.75,
                      verbose: bool = False
                      ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Resistant Fit Normalization to normalize the matrix counts.
    All the cells will be normalized to the centroid of the dataset.

    :param counts: Single Cell gene-count matrix. Rows are genes and columns are cells.
    :param p: the ratio of genes in the biological feature set.
    :param verbose: verbose flag for debug
    :return:
    norm_counts: normalized dataset (gene by cell)
    slope_list: the estimated slope for each cell
    intercept_list: the estimated intercept for each cell
    """
    if not isinstance(counts, np.ndarray):
        logger.warning(f"counter is not ndarray but {type(counts)}")
        if isinstance(counts, csr_matrix):
            logger.warning(f"counter is converted from csr_matrix to ndarray")
            counts = counts.toarray()
        else:
            raise NameError("counts must be either numpy.matrix or scipy.sparse.csr_matrix")

    counts = counts.astype('float32')  # float 32 should be enough based on scanpy normalization
    n_gene, n_cell = counts.shape
    logger.info(f'n_gene: {n_gene}')
    logger.info(f'n_cell: {n_cell}')

    norm_counts = counts.copy()
    slope_list = np.zeros(n_cell)
    intercept_list = np.zeros(n_cell)

    epoch = 1
    for j in range(epoch):
        logger.info(f'epoch: {j}')

        t_c = _get_centroid(norm_counts)
        logger.info(f'Seq-depth of centroid cell: {np.sum(t_c)}')
        for i in tqdm(range(n_cell)):
            res = _resistant_fit_linear(norm_counts[:, i], t_c, p=p, verbose=verbose)

            norm_counts[:, i] = res[0]
            slope_list[i] = res[1]
            intercept_list[i] = res[2]

    return norm_counts, slope_list, intercept_list
