import numpy as np


def geo_mean(iterable):
    a = np.array(iterable)
    # print('geo_mean:', len(a))
    return a.prod()**(1.0/len(a))


def median_of_ratio(x):
    # print()
    # x = x[np.sum(x, axis=1) >= 10, :]
    # print(x)
    x_g_mean = np.apply_along_axis(geo_mean, 1, x)
    # print('x_g_mean:\n', x_g_mean)
    # print('x_g_mean.shape:', x_g_mean.shape)
    # print('x_g_mean.shape:', x_g_mean[:, None].shape)
    x_divide = (x.T / x_g_mean).T
    # print('x_divide.shape: ', x_divide.shape)
    # print('x_divide:\n', x_divide)
    x_median_ratio = np.nanmedian(x_divide, axis=0)
    # print(x_median_ratio)
    # print(x_median_ratio.shape)
    res = x / x_median_ratio
    return res