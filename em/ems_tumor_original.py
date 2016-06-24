# -*- coding: UTF-8 -*-
"""Module containing the iterative em implementation.

"""
__version__ = '0.2'
__author__ = 'Esther Alberts'

import numpy as np

from utils import normalisation as norm
from utils import distributions as distr

import ems_initializer as ei

########################################################################


class IteraterEM(object):
    """Class to segment a brain image in tumor tissue and in brain tissues
    (white matter, grey matter and cerebrospinal fluid).

    Results are written to files in a subdirectory
    of the first input image (as supplied to InputData()).

    Attributes
    ----------
    invalid : boolean
        If true, the initializer is invalid and the em cannot be run.
    r : ems_read_input.ReaderEM
        Reader containing information about input data and run parameters.
    p : ems_initializer.Param
        Parameter instance.
    w : ems_write_output.WriterEM
        Writer instance.
    gauss_tissue_maps : 2D array
        Probabilistic tissue segmentations per gaussian
        (nb_gaussians, nr_brain_voxels).
    tumor_atlas : 2D array
        Specifies for every brain voxel how likeli it is to be tumor
        (first axis) or how likeli it is not to be tumor (second axis).
    tumor_maps : 2D array
        Probabilistic tumor segmentation for brain voxels
        (nr_channels, nr_brain_voxels).
    means : 2D array
        means per gaussian used to model tissue intensity distribution
        per channel (nr_gaussians, nr_channels).
    variances : ndarray
        variances per gaussian used to model tissue intensity distribution
        per channel (nr_gaussians, nr_channels).
    tumor_means
        means per gaussian used to model tumor intensity distribution
        per channel (nr_gaussians, nr_channels).
    tumor_vars
        variances per gaussian used to model tumor intensity distribution
        per channel (nr_gaussians, nr_channels).
    converged : boolean
        If true, algorithm has finished its iterations.
    iteration : int
        Iteration number.
    log_llh_arr : 1D array
        Log likelihood values per iteration.
    log_llh : float
        Most recent Log likelihood.
    rel_log_llh : float
        Relative increase in Log likelihood to previous iteration.

    """

    def __init__(self, em_initializer):
        """
        Initialises all variables needed to call self.run().

        Parameters
        ----------
        em_initializer : ems_initializer.Intializer
            Initialisation instance.

        """
        if not isinstance(em_initializer, ei.Initializer):
            raise ValueError()

        if not em_initializer.invalid:

            self.invalid = False

            self.r = em_initializer.r
            self.p = em_initializer.p
            self.w = em_initializer.w

            self.gauss_tissue_maps = self.r.gauss_tissue_maps_init
            self.tumor_maps = self.r.tumor_maps_init
            self.tumor_atlas = self.r.tumor_atlas_init

            self._em_initialisation()

        else:

            print 'Invalid initializer, this object is obsolete'

            self.invalid = True

    def _em_initialisation(self):
        """Called by self.__init__() to initialise variables."""
        # Initialisations for intensity guassians.
        self.means = np.zeros(shape=(self.r.nr_gaussians,
                                     self.r.nr_channels))
        self.variances = np.zeros(shape=(self.r.nr_gaussians,
                                         self.r.nr_channels,
                                         self.r.nr_channels))
        self.tumor_means = np.zeros(shape=(self.r.nr_channels,))
        self.tumor_vars = np.zeros(shape=(self.r.nr_channels,))

        # Initialisations for iterations.
        self.converged = 0
        self.iteration = 0

        # Log-likelihood initialisations.
        self.log_llh_arr = np.zeros((self.p.max_iteration,))
        self.log_llh = 1 / np.float(ei.EPS)
        self.rel_log_llh = 1

    ####################################################################

    def run(self):
        """Start em iterations """

        self._print_start_message()

        # Start iterating.
        while self.converged == False:

            self.iteration = self.iteration + 1
            self.old_log_llh = self.log_llh

            self._print_start_iteration()

            ##########################
            # 1. Calculate normal means.
            self._update_tissue_means()
            # 2. Calculate tumor means.
            self._update_tumor_means()
            # 3. Calculate classification.
            self._evaluate_configurations()
            # 4. Calculating normal tissue classification
            #    and tumor classification.
            self._update_classification_maps()
            ###########################

            # Convergence criterion
            self.converged = (self.iteration == self.p.max_iteration)

        self.w.write_tissue_maps(self.tissue_maps)
        self.w.write_tumor_maps(self.tumor_maps, self.tumor_atlas)


    def _print_start_message(self):
        """Print a start message."""

        print('')
        print('Simple EMS has started')
        print(' - no biasfield correction')
        print(' - no MRF spatial regularisator')
        print('**Input data is supposed to be preprocessed:**')
        print('Images, masks and atlas should...')
        print(' - be in the same reference space')
        print(' - have the same pixel resolution')
        print(' - have the same dimensions')
        print('*No tumor will be allowed in CSF voxels*')
        print('')
        print('FLAT PRIOR: ' + str(self.p.flat_prior))

    def _print_start_iteration(self):
        """Print a start message before iterating."""

        print('---')
        print('Iteration ' + str(self.iteration))
        print('')

    def _update_tissue_means(self):
        """Function calculating gaussian parameters for healthy tissue.

        Call self.gauss_tissue_maps, self.tumor_maps, self.r.data.
        Modify self.means, self.variances.
        
        """
        print('Estimating gaussian parameters for healthy tissue')
        print('')
        for gauss in range(self.r.nr_gaussians):

            for channel in range(self.r.nr_channels):
                weights = self.gauss_tissue_maps[gauss, :] * \
                            (1 - self.tumor_maps[channel, :])

                if np.sum(weights) == 0:
                    err = 'No voxels remain to update gaussian for ' + \
                          'this tissue: tissue weights are all zero'
                    raise ZeroDivisionError(err)

                self.means[gauss, channel] = \
                    distr.get_mean(self.r.data[channel, :], weights)
                self.variances[gauss, channel, channel] = \
                    distr.variance(self.r.data[channel, :],
                                   weights,
                                   self.means[gauss, channel])

        # Split clusters within the same tissue class:
        # they are initialised identically and will remain identical
        # unless we push them apart artificially
        if self.iteration <= 3:
            print('Splitting clusters of normal tissue classes')
            print(' - old means:')
            print(str(self.means))
            self.means = split_class_means(self.means, self.r.lkp)

        print(' - normal means of iteration ' + str(self.iteration) + ': ')
        print(str(self.means))
        print('')

    def _update_tumor_means(self):
        """Function calculating gaussian parameters for healthy tissue.

        Call self.tumor_maps.
        Modify self.tumor_weights, self.tumor_means, self.tumor_vars.
        
        """

        print('Estimating gaussian parameters for tumor tissue')
        print('')
        new_tumor_weights = np.zeros([self.r.nr_channels,
                                      self.r.nr_brain_voxels])
        for channel in range(self.r.nr_channels):

            weights = self.tumor_maps[channel, :]

            weights = self._select_tumor_intensities(weights, channel)

            if np.sum(weights) == 0:
                err = 'No tumor pixels present to update tumor gaussian'
                err += '\nIs tumor present in the input images?'
                err += '\nIf so, try loosening the tumor threshold'
                err += '\n ( Check `flat_prior` in code )'
                print err
                if self.iteration == 1:
                    new_tumor_weights[channel, :] = 0
                    # Remove all configurations with tumor in this channel
                    channel_column = self.r.tumor_transitions[:, channel]
                    tumor_rows = np.nonzero(channel_column == ei.TUMOR)[0]
                    for row in tumor_rows:
                        np.delete(self.r.tumor_transitions, row, 0)
                        np.delete(self.r.gaussians, row, 0)
                        np.delete(self.r.classes, row, 0)
                    self.r.nr_tissue_conf = \
                        self.r.tumor_transitions.shape[0]
                else:
                    new_tumor_weights[channel, :] = \
                        np.copy(self.tumor_weights[channel, :])
            else:
                new_tumor_weights[channel, :] = weights

            self.tumor_means[channel] = \
                distr.get_mean(self.r.data[channel, :],
                               new_tumor_weights[channel, :])
            self.tumor_vars[channel] = \
                distr.variance(self.r.data[channel, :],
                               new_tumor_weights[channel, :],
                               self.tumor_means[channel])

        self.tumor_weights = new_tumor_weights

        print('Tumor means of iteration ' + str(self.iteration) + ': ')
        print(str(self.tumor_means))
        print('')

    def _select_tumor_intensities(self, weights, channel):
        """Select voxels with intensities likeli to be tumor based on 
        the current gaussians modeling the healthy tissues.

        Calls self.means, self.variances, self.r.lkp, self.r.data.
        
        """
        if np.sum(weights) != 0:
            channel_means = self.means[:, channel]
            mean_gm = np.average(channel_means[self.r.lkp == ei.GM])
            mean_wm = np.average(channel_means[self.r.lkp == ei.WM])
            # distance between white and gray matter means
            wm_gm_diff = np.absolute(mean_gm - mean_wm)

            for gauss in np.nditer(np.where(self.r.lkp == ei.WM)):
                wm_prob = distr.gaussian(\
                            self.r.data[channel, :],
                            self.means[gauss, channel],
                            self.variances[gauss, channel, channel],
                            False)
                wm_prob = np.amax(\
                    [wm_prob, np.zeros((self.r.nr_brain_voxels,))], axis=0)

            not_wm_prob = 1 - wm_prob
            not_wm_prob = norm.window(not_wm_prob, 0, 1)
            weights = weights * not_wm_prob

            # Throw out pixels which are far from hyper- or
            # hypo-intensive depending on the modality type
            if self.r.hyper[channel] == 1:
                mple = 1
                if self.p.t1c_extra_hyper and self.r.is_t1c_channel(channel):
                    mple = 3
                weights[self.r.data[channel, :] <= \
                    (self.means[ei.WM, channel] + (mple * wm_gm_diff))] = 0
            else:
                weights[self.r.data[channel, :] >= \
                    (self.means[ei.WM, channel] - wm_gm_diff)] = 0

        return weights

    def _evaluate_configurations(self):
        """Evaluate configurations incorporating prior and
        current parameters.

        Call self.r.tumor_transitions, self.r.atlas_prior,
        self.r.guassians, self.r.tumor_transitions,
        self.r.data, self.tumor_atlas,
        self.means, self.variances, self.tumor_means, self.tumor_vars.

        Modify self.tissue_tumor_maps, self.log_llh,
        self.log_llh_arr, self.rel_log_llh.

        """
        print('Estimating tissue/tumor segmentations')
        print('')
        likeli_product = np.zeros(shape=(self.r.nr_brain_voxels,))
        tissue_tumor_maps = np.zeros(shape=(self.r.nr_tissue_conf,
                                            self.r.nr_brain_voxels))
        for conf in range(self.r.nr_tissue_conf):

            gauss = self.r.gaussians[conf]
            tumor_array = self.r.tumor_transitions[conf, :]

            likeli_product[:] = self.r.atlas_prior[gauss, :]

            if self.p.flat_prior == 1 and self.iteration == 1:
                self._update_tumor_early()

            for channel in range(self.r.nr_channels):
                # If no tumor in channel.
                if tumor_array[channel] == ei.NO_TUMOR:
                    mean = self.means[gauss, channel]
                    var = self.variances[gauss, channel, channel] + ei.EPS

                    likeli_product = likeli_product * \
                        distr.gaussian(self.r.data[channel, :],
                                       mean,
                                       var,
                                       self.p.gaussian_norm)

                    likeli_product = likeli_product * \
                        self.tumor_atlas[ei.NO_TUMOR, :]
                # If tumor in channel.
                if tumor_array[channel] == ei.TUMOR:
                    mean = self.tumor_means[channel]
                    var = np.amax(self.tumor_vars[channel], ei.EPS)

                    # artificially ensure that more intense voxels
                    # in the hyperintense channels are not less likeli
                    # to be tumor (caveat of evaluating a nonsymmetrical
                    # distribution with a gaussian)
                    adapted_data = np.copy(self.r.data[channel, :])
                    if self.r.hyper[channel] and \
                            (self.iteration >= self.p.max_iteration - 2):
                        data_tmp = adapted_data
                        data_tmp = data_tmp - np.sqrt(var)
                        data_tmp[data_tmp < 0] = 0
                        data_tmp[data_tmp > mean] = mean
                        adapted_data = data_tmp

                    likeli_product = likeli_product * distr.gaussian(\
                        adapted_data, mean, var, self.p.gaussian_norm)

                    likeli_product = likeli_product * \
                        self.tumor_atlas[ei.TUMOR, :]

            tissue_tumor_maps[conf, :] = np.fmax(likeli_product, 0)

        # Calculate likelihood and normalise segmentation.
        likelihood = np.fmax(np.sum(tissue_tumor_maps, axis=0),
                             ei.TINY)
        self.tissue_tumor_maps = tissue_tumor_maps / \
            np.tile(likelihood, (self.r.nr_tissue_conf, 1))
        self.log_llh = np.sum(np.log(likelihood))
        self.log_llh_arr[self.iteration - 1] = self.log_llh
        self.rel_log_llh = np.absolute(\
            (self.old_log_llh - self.log_llh) / self.log_llh)

        print('=> logLikelihood = ' + str(self.log_llh))
        print('=> relative change in loglikelihood = ' + \
              str(self.rel_log_llh))
        print('')

    def _update_tumor_early(self):
        """Update the flat tumor initialisation already
        with the current gaussian parameters.

        Modify self.tumor_atlas, self.tumor_maps.
        """
        for channel in range(self.r.nr_channels):
            mean = self.tumor_means[channel]
            var = np.amax(self.tumor_vars[channel], ei.EPS)
            self.tumor_maps[channel] = distr.gaussian(\
                self.r.data[channel, :], mean, var, self.p.gaussian_norm)
        self.tumor_atlas[ei.TUMOR, :] = \
            np.sum(self.tumor_maps, axis=0) / self.r.nr_channels
        self.tumor_atlas[ei.NO_TUMOR, :] = \
            1 - self.tumor_atlas[ei.TUMOR, :]

    def _update_classification_maps(self):
        """Get the new tissue and tumor segmentations
        or this iteration.

        Calls self.tissue_tumor_maps, self.r.tumor_transitions
        self.r.tumor_mask, self.r.gaussians, self.r.classes.

        Modifies self.tumor_atlas, self.tumor_maps,
        self.gauss_tissue_maps, self.tissue_maps.

        """
        print('Calculating normal tissue and tumor classification')
        print('')
        tissue_maps = np.zeros(shape=(self.r.nr_classes,
                                      self.r.nr_brain_voxels))
        gauss_tissue_maps = np.zeros(shape=(self.r.nr_gaussians,
                                            self.r.nr_brain_voxels))
        tumor_maps = np.zeros(shape=(self.r.nr_channels,
                                     self.r.nr_brain_voxels))

        for conf in range(self.r.nr_tissue_conf):
            configuration_prob = self.tissue_tumor_maps[conf, :]
            tumor_array = self.r.tumor_transitions[conf, :]

            tumor_channels = np.where(tumor_array == ei.TUMOR)
            tumor_maps[tumor_channels, :] += \
                    configuration_prob * self.r.tumor_mask

            class_nr = self.r.classes[conf]
            gauss = self.r.gaussians[conf]
            tissue_maps[class_nr, :] += configuration_prob
            gauss_tissue_maps[gauss, :] += configuration_prob

        tissue_maps = norm.window(tissue_maps, 0, 1)
        gauss_tissue_maps = norm.window(gauss_tissue_maps, 0, 1)
        tumor_maps = norm.window(tumor_maps, 0, 1)

        # Update tumor atlas.
        self.tumor_atlas[ei.TUMOR, :] = \
            np.sum(tumor_maps, axis=0) / self.r.nr_channels
        self.tumor_atlas[ei.NO_TUMOR, :] = \
            1 - self.tumor_atlas[ei.TUMOR, :]

        self.tumor_maps = tumor_maps
        self.gauss_tissue_maps = gauss_tissue_maps
        self.tissue_maps = tissue_maps


########################################################################


def split_class_means(means, class_per_gaussian):
    """Split means within the same channel and class
    when they have the same mean value.

    Parameters
    ----------
    means: 2D array
        Means per gaussian (1st axis) and per channel (2nd axis).
    class_per_gaussian: 1D array
        Indices for each gaussian to which tissue class it belongs.

    """

    changes = [1.005, 0.995, 1.007, 0.993, 1.01, 0.99]
    nr_channels = means.shape[1]
    nr_priors = len(np.unique(class_per_gaussian))
    for channel in range(nr_channels):
        channel_means = means[:, channel]
        for class_nr in range(nr_priors):
            class_means = channel_means[class_per_gaussian == class_nr]
            unique_means = np.unique(class_means)
            nonuniques = (len(class_means) != len(unique_means))
            while nonuniques:
                for mean in unique_means:
                    if np.count_nonzero(class_means == mean) > 1:
                        ind = (class_means == mean)
                        class_means[ind] = \
                            [x * mean
                             for x in changes[0:np.count_nonzero(ind)]]
                unique_means = np.unique(class_means)
                nonuniques = (len(class_means) != len(unique_means))
            channel_means[class_per_gaussian == class_nr] = class_means
        means[:, channel] = channel_means

    return means


########################################################################
