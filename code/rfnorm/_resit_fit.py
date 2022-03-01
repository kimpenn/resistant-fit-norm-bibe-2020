import logging.config
from typing import Tuple

import numpy as np
from scipy import stats
from sklearn import linear_model

logging.config.fileConfig(fname='logger.conf', disable_existing_loggers=False)
logger = logging.getLogger(__name__)


def _preprocess(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Return the index of genes which are valid for regression model.
    The following outlier genes will be excluded since they have no contribution to infer the model.
    - If the gene count is zero in at least one cell.
    - If the gene has extremely high expression level which will break the regression model.

    By default, extremely high expression means: rank < 0.01 among non-zero genes.

    :param x: the gene expression of the cell x, reference cell
    :param y: the gene expression of the cell y
    :return: boolean mask of the genes. True for valid genes and False for the outliers.
    """
    # TODO:: Fact check the outlier rate
    outlier_rate = 0.01
    select_mask = np.logical_or(x != 0, y != 0).reshape(-1)
    n_select = np.sum(select_mask)
    clip_idx = np.flip(np.argsort(x + y))[0:int(n_select * outlier_rate)]
    select_mask[clip_idx] = False

    return select_mask


def _resistant_fit_linear(
        source: np.ndarray,
        target: np.ndarray,
        p: np.float32 = 0.75,
        verbose: bool = False,
        init_step: str = "ransac"
) -> Tuple[np.ndarray,
           np.float32,
           np.float32]:
    """
    Use resistant fit to normalize cell x (source) to the reference cell y (target).

    :param source: the gene expression of the cell x
    :param target: the gene expression of the cell y, reference cell
    :param p: the size of biological feature set
    :param verbose: verbose flag for debug
    :return:
        y_regression: the normalized 1d array
        slope: the final slope from EM Regression
        intercept: the final intercept from EM Regression
    """
    ########################################
    # Select valid genes for regression
    ########################################
    iter_limit = 20
    np.seterr(all='raise')
    select_mask = _preprocess(source, target)
    n_select = np.sum(select_mask)

    x_select = source[select_mask].copy()  # Note that len(x_select) <= source
    y_select = target[select_mask].copy()

    ############################################################
    # Init EM step: robust regression on all genes
    ############################################################
    # Robust regression on the whole dataset to ignore outliers
    if init_step == "ransac":
        ransac = linear_model.RANSACRegressor(random_state=42)
        ransac.fit(x_select.copy().reshape(-1, 1), y_select.copy().reshape(-1, 1))
        slope, intercept = float(ransac.estimator_.coef_), float(ransac.estimator_.intercept_)
    elif init_step == "siegel":
        slope, intercept = stats.siegelslopes(y_select, x_select)
    elif init_step == "theil":
        slope, intercept, _, _ = stats.theilslopes(y_select, x_select)
    else:
        raise NameError("init_step must be chosen from list ['ransac', 'siegel', 'theil']")

    y_regression = np.asarray([slope * x_iter + intercept for x_iter in x_select])
    square_list = np.square(y_select - y_regression)
    square_list_index_sort = np.argsort(square_list)
    sub_index = square_list_index_sort[0:int(n_select * p) + 1]

    # Set Biological Feature Set (BFS)
    x_bfs = x_select[sub_index]
    y_bfs = y_select[sub_index]

    if verbose:
        logger.info(f'[Init EM step] slope:{slope},  intercept:{intercept}')

    ############################################################
    # Resistant Fit Regression on BFS
    ############################################################
    loss_pre = np.Inf
    for i in range(iter_limit):
        # E step: Linear regression on BFS
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_bfs, y_bfs)

        y_regression = np.asarray([slope * x_iter + intercept for x_iter in x_select])
        square_list = np.square(y_select - y_regression)
        square_list_index_sort = np.argsort(square_list)
        sub_index = square_list_index_sort[1:int(n_select * p) + 1]

        loss = np.sum(square_list[sub_index])

        delta_loss = abs(loss_pre - loss)
        if delta_loss < 0.00001:
            if verbose:
                logger.info('convergence')
            break
        else:
            loss_pre = loss

        # M step: Update x and y for next iteration of linear regression
        x_bfs = x_select[sub_index]
        y_bfs = y_select[sub_index]

        if i == iter_limit:
            if not verbose:
                logger.info('[Resist Fit] Reach iteration limit')

    ############################################################
    # Normalize cell y based on regression model
    ############################################################
    y_regression = slope * source + intercept
    # If x is zero then clip y to 0
    # TODO: x is unlikely to be zero.
    y_regression[source == 0] = 0

    if verbose:
        logger.info(f'[Resist Fit] slope: {slope}, intercept: {intercept}')
        logger.info(f'depth(y_select): {np.sum(y_select)}, \n'
                    f'depth(x_select): {np.sum(x_select)}, \n'
                    f'depth(x): {np.sum(source)},\n'
                    f'depth(norm): {np.sum(y_regression)},\n'
                    f'depth(y): {np.sum(target)}\n'
                    f'y_select/x_select: {np.sum(y_select) / np.sum(x_select)}')

    return np.float32(y_regression), np.float32(slope), np.float32(intercept)


# def _resistant_fit_ransac(
#         source: np.ndarray,
#         target: np.ndarray,
#         p: np.float32 = 0.75,
#         verbose: bool = False
# ) -> np.ndarray:
#     ########################################
#     # Select valid genes for regression
#     ########################################
#     np.seterr(all='raise')
#     select_mask = _preprocess(source, target)
#     x_select = source[select_mask].reshape(-1, 1)
#     y_select = target[select_mask].reshape(-1, 1)
#     x_res = source.copy()
#
#     ########################################
#     # Robust Regression: RANSAC
#     ########################################
#     ransac = linear_model.RANSACRegressor(random_state=42)
#     ransac.fit(x_select, y_select)
#     x_res[select_mask] = ransac.predict(x_select)[:, 0]
#
#     return x_res
