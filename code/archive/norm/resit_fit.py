import numpy as np
from scipy import stats
import statsmodels.api as sm
from pandas import Series, DataFrame


def poly_resist_fit(x, y):
    print('Running poly resist fit...')
    n_gene = x.shape[0]
    print('number of genes: ', n_gene)
    x_reg = x
    y_reg = y
    loss_pre = float('Inf')

    # Set the initial points of slope and intercept
    for i in range(3000):
        z = np.polyfit(x_reg, y_reg, 3)
        p_fit = np.poly1d(z)
        y_hat = p_fit(x)
        square_list = np.square(y - y_hat)
        square_list_index_sort = np.argsort(square_list)
        sub_index = square_list_index_sort[1:int(n_gene / 2) + 1]

        loss = np.sum(square_list[sub_index])

        delta_loss = abs(loss_pre - loss)
        if delta_loss < 0.00001:
            print('convergence')
            break
        else:
            loss_pre = loss
        # Update x and y for next iteration of linear regression
        x_reg = x[sub_index]
        y_reg = y[sub_index]
        print(i, loss, delta_loss)

    return p_fit, y_hat


def linear_resist_fit(x, y):
    np.seterr(all='raise')
    # print('Running linear resist fit...')
    n_gene = x.shape[0]
    # print('number of genes: ', n_gene)
    x_reg = x
    y_reg = y
    loss_pre = float('Inf')
    # Set the initial points of slope and intercept
    for i in range(3000):
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_reg, y_reg)
        abline_values = np.asarray([slope * x_iter + intercept for x_iter in x])
        square_list = np.square(y - abline_values)
        square_list_index_sort = np.argsort(square_list)
        # sub_index = square_list_index_sort[1:int(n_gene * 0.3) + 1]
        sub_index = square_list_index_sort[1:int(n_gene / 2) + 1]

        loss = np.sum(square_list[sub_index])

        delta_loss = abs(loss_pre - loss)
        if delta_loss < 0.00001:
            print('convergence')
            break
        else:
            loss_pre = loss
        # Update x and y for next iteration of linear regression
        x_reg = x[sub_index]
        y_reg = y[sub_index]
        # print(i, loss, delta_loss)

    return abline_values


def linear_resist_fit_p(x, y, p=0.5):
    np.seterr(all='raise')
    # print('Running linear resist fit...')
    n_gene = x.shape[0]
    # print('number of genes: ', n_gene)
    x_reg = x
    y_reg = y
    loss_pre = float('Inf')
    # Set the initial points of slope and intercept
    for i in range(3000):
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_reg, y_reg)
        abline_values = np.asarray([slope * x_iter + intercept for x_iter in x])
        square_list = np.square(y - abline_values)
        square_list_index_sort = np.argsort(square_list)
        # sub_index = square_list_index_sort[1:int(n_gene * 0.3) + 1]
        sub_index = square_list_index_sort[1:int(n_gene * p) + 1]

        loss = np.sum(square_list[sub_index])

        delta_loss = abs(loss_pre - loss)
        if delta_loss < 0.00001:
            # print('convergence')
            break
        else:
            loss_pre = loss
        # Update x and y for next iteration of linear regression
        x_reg = x[sub_index]
        y_reg = y[sub_index]
        # print(i, loss, delta_loss)

    return abline_values


def linear_resist_fit_robust(x, y, p=0.5):
    np.seterr(all='raise')
    # print('Running linear resist fit...')
    n_gene = x.shape[0]
    # print('number of genes: ', n_gene)
    x_reg = x
    y_reg = y
    loss_pre = float('Inf')
    # Set the initial points of slope and intercept
    for i in range(3000):
        slope, intercept = stats.siegelslopes(x_reg, y_reg)
        abline_values = np.asarray([slope * x_iter + intercept for x_iter in x])
        square_list = np.square(y - abline_values)
        square_list_index_sort = np.argsort(square_list)
        # sub_index = square_list_index_sort[1:int(n_gene * 0.3) + 1]
        sub_index = square_list_index_sort[0: int(n_gene * p) + 1]

        loss = np.sum(square_list[sub_index])

        delta_loss = abs(loss_pre - loss)
        if delta_loss < 0.00001:
            print('convergence')
            break
        else:
            loss_pre = loss
        # Update x and y for next iteration of linear regression
        x_reg = x[sub_index]
        y_reg = y[sub_index]
        print(i, loss, delta_loss)

    return abline_values


def select_zero_genes(x, y):
    return np.multiply(x, y) == 0


def preprocess(x, y, p_default=0.5):
    """
    Estimate the size of biological set based on the input data.
    Assume data is gene(row) x cell(col) matrix.

    By default, 50% of genes is a reasonable choice of size.
    But if the number of zeros in one column is more than 50%,
    program will cast warning and set p = 0.7 * min(non-zeros).

    :param y:
    :param x:
    :param p_default:
    :return: the estimated size p
    """

    def select_genes(x, y):
        outlier_rate = 0.01
        select_list = np.multiply(x, y) > 0
        n_select = np.sum(select_list)
        clip_idx = np.flip(np.argsort(x+y))[0:int(n_select * outlier_rate)]
        select_list[clip_idx] = False

        return select_list
    
    n_gene = y.shape[0]
    select_list = select_genes(x, y)

    p = 0.9
    # p = 0.7
    # if n_select_list < n_gene * p_default:
    #     # Shrink of size of biological feature set
    #     print("warning: Shrink of size of biological feature set")
    #     p = 0.7
    #     print(f'p: {p}')
    # else:
    #     p = p_default
    return select_list, p


def linear_resist_fit_robust_mix(x, y, p_default=0.5, verbose=False):
    iter_limit = 20
    np.seterr(all='raise')
    # print('Running linear resist fit...')
    n_gene = x.shape[0]
    # print('number of genes: ', n_gene)
    select_list, p = preprocess(x, y)
    # print(f'p: {p}')
    n_select = np.sum(select_list)

    x_select = x[select_list].copy()
    y_select = y[select_list].copy()

    x_reg = x_select
    y_reg = y_select
    loss_pre = float('Inf')

    # Robust regression on the whole dataset to ignore outliers
    slope, intercept = stats.siegelslopes(y_reg, x_reg)
    abline_values = np.asarray([slope * x_iter + intercept for x_iter in x_select])
    square_list = np.square(y_select - abline_values)
    square_list_index_sort = np.argsort(square_list)
    sub_index = square_list_index_sort[0:int(n_select * p) + 1]
    x_reg = x_select[sub_index]
    y_reg = y_select[sub_index]
    if verbose:
        print(f'[siegelslopes] slope:{slope},  intercept:{intercept}')
    # Set the initial points of slope and intercept
    for i in range(iter_limit):
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_reg, y_reg)
        # slope, intercept = stats.siegelslopes(y_reg, x_reg)

        abline_values = np.asarray([slope * x_iter + intercept for x_iter in x_select])
        square_list = np.square(y_select - abline_values)
        square_list_index_sort = np.argsort(square_list)
        sub_index = square_list_index_sort[1:int(n_select * p) + 1]

        loss = np.sum(square_list[sub_index])

        delta_loss = abs(loss_pre - loss)
        if delta_loss < 0.00001:
            if verbose:
                print('convergence')
            break
        else:
            loss_pre = loss
        # Update x and y for next iteration of linear regression
        x_reg = x_select[sub_index]
        y_reg = y_select[sub_index]

        if i == iter_limit:
            if not verbose:
                print('[Resist Fit] Reach iteration limit')

    # abline_values = np.asarray([slope * x_iter + intercept for x_iter in x])
    abline_values = slope * x + intercept
    # abline_values = slope * x
    # print(f'# of 0s in abline: {np.sum(abline_values == 0)}')
    # abline_values[select_zero_genes(x, y)] = 0
    abline_values[x == 0] = 0
    # print(f'# of 0s in new abline: {np.sum(abline_values == 0)}')
    if verbose:
        print(f'[Resist Fit] slope: {slope}, intercept: {intercept}')
        print(f'depth(y_select): {np.sum(y_select)}, \n'
          f'depth(x_select): {np.sum(x_select)}, \n'
          f'depth(x): {np.sum(x)},\n'
          f'depth(norm): {np.sum(abline_values)},\n'
          f'depth(y): {np.sum(y)}\n'
          f'y_select/x_select: {np.sum(y_select) / np.sum(x_select)}')
    return abline_values, slope, intercept
