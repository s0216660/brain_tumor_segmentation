# -*- coding: UTF-8 -*-
"""Module including basic functions for data normalisation.

"""
__version__ = '0.2'
__author__ = 'Esther Alberts'

import numpy as np

########################################################################

def window(data, lower=0, upper=1):
    """Window `data` with `lower` and `upper` thresholds."""

    data[data < lower] = lower
    data[data > upper] = upper
    return data

def rescale(data, lower=0, upper=1):
    """Scale `data` between `lower` and `upper`."""

    norm = float(np.amax(data) - np.amin(data))
    if norm != 0:
        new_data = ((upper - lower) * \
                    (data - np.amin(data)) / norm) + lower
    else:
        new_data = (data - np.amin(data)) + lower

    return new_data

def whiten(data):
    """Whiten `data` to obtain zero means and unit std."""

    norm = float(np.std(data))
    if norm != 0:
        new_data = (data - np.mean(data)) / norm
    else:
        new_data = data - np.mean(data)

    return new_data
