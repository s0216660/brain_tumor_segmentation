# -*- coding: UTF-8 -*-
"""Module allowing to generate and evaluate probabilistic distributions

"""
__version__ = '0.2'
__author__ = 'Esther Alberts'

import numpy as np

########################################################################

EPS = np.finfo(np.double).eps

########################################################################

def logistic(data, mean, var, normalized=True):
    """Logistic evaluation of `data`, `normalized` if wished,
    given `mean` and `var`.

    """

    if var == 0:
        err = 'Gaussian evaluation fails: zero mean'
        raise ZeroDivisionError(err)

    neg_exp = np.exp(-(data - mean) / np.float(np.sqrt(var)))
    density = 1 / (1 + neg_exp)
    if normalized:
        norm = (np.sqrt(2 * np.pi * var))
        density = density / norm
    return density

def gaussian(data, mean, var, normalized=True):
    """Get density of `data` under a (`normalized`) gaussian with
    parameters `mean` and `var`.

    """

    if var == 0:
        err = 'Gaussian evaluation fails: zero variance'
        raise ZeroDivisionError(err)

    mahalanobis_sq = (data - mean) ** 2 / np.float(var)
    density = np.exp(-mahalanobis_sq / 2)
    if normalized:
        density = density / (np.sqrt(2 * np.pi * var))
    return density

def frechet(data, mean, std, scale):
    """Evaluate `data` under a frechet distribution with
    parameters `mean`, `std` and `scale`.

    """

    if std == 0:
        err = 'Frechet evaluation fails: zero std'
        raise ZeroDivisionError(err)

    density = np.zeros_like(data)
    distance = (data[data >= mean] - mean) / std
    density[data >= mean] = (scale / std) * \
                          (distance ** (-scale - 1)) * \
                          np.exp(-(distance ** (-scale)))

    return density

def gaussian_multivar(data, means, covars, normalized):
    """Calculate overall [1 x s] `density` of a [n x s] `data` under a
    n-multivariate gaussian with an array of [n x 1] `means` and a
    [n x n] covariance matrix `covars`.

    """

    dim_data = len(data.shape)

    # In case the gaussian is univariate (n=1)
    if dim_data == 1:
        return gaussian(data, means, covars, normalized)

    else:
        nr_samples = data.shape[1]
        nr_features = data.shape[0]

        data_diff = data - np.outer(means, np.ones(shape=(nr_samples)))
        det_covar = np.linalg.det(covars)

        if det_covar == 0:
            err = 'Zero determinant of the mean matrix'
            raise ValueError(err)

        mahalanobis_sq = np.sum(np.dot(np.transpose(data_diff),
                                       np.linalg.inv(covars))
                                * np.transpose(data_diff),
                                axis=1)

        if np.any(mahalanobis_sq < 0):
            err = 'Negative mahalanobis distances'
            raise ValueError(err)
        density = np.exp(-mahalanobis_sq / 2)

        if normalized:
            norm = np.sqrt(((2 * np.pi) ** nr_features) * det_covar)
            density = density / norm

    return density

def get_mean(data, weights):
    """Calculate the mean of the [n x s] `data` given a [n x s] set of
    `weights`.

    """

    if len(weights.shape) == 1 and len(data.shape) == 1:
        sum_axis = 0
    else:
        sum_axis = 1

    norm = np.sum(weights, axis=sum_axis)
    mean = np.zeros_like(norm)

    if isinstance(norm, np.float) or isinstance(norm, float):
        if norm == 0:
            return np.NaN
        else:
            return np.sum((weights * data), axis=sum_axis) / norm
    else:
        mean[norm == 0] = np.NaN
        mean[norm != 0] = np.sum((weights * data),
                                 axis=sum_axis)[norm != 0] \
                        / norm[norm != 0]

    return mean

def variance(data, weights, mean):
    """Calculate variance of a [s,] dataset together with a [s,] set
    of weights and a mean.

    Parameters
    ----------
    data: 1D array [s,]
        dataset of s samples and n channels
    weights: 1D array [s,]
        weights corresponding to the data
    means: np.float
         mean corresponding to the data

    Returns
    -------
    - var: np.float
        variance of the data given the sets of weights and
        the mean of the data.

    """

    if np.count_nonzero(weights < 0) is not 0:
        warning = '**Negative weights present for the gaussian update**'
        print warning

    data_diff = data - mean
    norm = np.sum(weights)
    if norm == 0:
        var = 0
    else:
        var = np.sum((data_diff ** 2) * weights) / norm

    return np.squeeze(var)

def covariance(data, weights, means, verbose=True):
    """Calculate mean of a [n x s] dataset together with a [n x s]
    set of weights and [n x 1] means.

    Parameters
    ----------
    data: 2D array [n x s]
        dataset of s samples and n channels
    weights: 2D array [n x s]
        weights corresponding to the datasets in data
    means: 1D array [1 x n]
        array of means corresponding to the datasets in data

    Returns
    -------
    var: mean matrix of the datasets given the sets of weights and the
    array of means corresponding to these datasets.

    """

    if np.count_nonzero(weights < 0) is not 0:
        warning = '**Negative weights are not allowed**'
        print warning
    nonsymm = False

    if len(weights.shape) == 1 and len(data.shape) == 1:
        if verbose:
            print 'You were looking for the variance, not for the mean'
        return variance(data, weights, means)

    elif len(weights.shape) != len(data.shape):
        err = 'Data and weights are not of compatible dimension'
        raise AttributeError(err)

    else:
        data_diff = data - np.outer(means,
                                    np.ones(shape=(data.shape[1])))

        # check if weights overlap,
        # if not assume variables are independent
        norms = np.dot(np.sqrt(weights),
                       np.transpose(np.sqrt(weights)))
        norms[norms == 0] = 1

        var = np.dot(data_diff * np.sqrt(weights),
                     np.transpose(np.sqrt(weights) * data_diff)) \
              / norms

        if not np.all(var == np.transpose(var)):
            var = (var + np.transpose(var)) / 2
        if not np.all(var == np.transpose(var)):
            warning = '**Covariance matrix remains asymmetrical.\n' + \
                      'Ignore if no further problems arise**'
            print warning
            nonsymm = True

        if np.any(np.diag(var) < 0):
            err = 'negative variances on diagonal of mean matrix'
            raise ValueError(err)

        eig, eigv = np.linalg.eigh(var)
        if np.any(eig <= np.sqrt(EPS)):
            print('Covariance matrix is not pos semidef')
            print('Neg eigenvalues are corrected...')
            var = nearest_PSD_matrix(var, eig, eigv, verbose)

    return np.squeeze(var), nonsymm

def nearest_PSD_matrix(covars, eig, eigv, verbose=True):
    """Calculate the nearest pos semi-definite matrix of the given
    covars matrix with eigenvalues eig and eigenvectors eigv.

    """

    if verbose:
        print('Non-pos-semidef mean matrix: ' + str(covars))
        print('Eigenvalues ' + str(eig))

    Q = np.matrix(eigv)
    xdiag = np.matrix(np.diag(np.maximum(eig, np.sqrt(EPS))))
    covars = Q * xdiag * Q.T

    new_eig = np.linalg.eigh(covars)[0]
    if verbose:
        print('Corrected pos-semidef mean matrix: ' + str(covars))
        print('New eigenvalues ' + str(new_eig))

    if np.any(new_eig <= 0):
        err = 'After correction of covar matrix still negative or '\
              'zero eigenvalues!'
        raise ValueError(err)

    return covars
