# -*- coding: UTF-8 -*-
"""Module containing
- two classes containing parameters: Param() and TumorLocation()
- class for initialisation of em brain tumor segmentation

"""
__version__ = '0.2'
__author__ = 'Esther Alberts'

import numpy as np

from utils import mr_variables as mr

import ems_read_input as ri
import ems_write_output as wo

########################################################################

GM = 0
WM = 1
CSF = 2
TISSUES = ['gm', 'wm', 'csf']

EPS = np.finfo(np.double).eps
TINY = np.finfo(np.double).tiny

TUMOR = 1
NO_TUMOR = 0

########################################################################


class Param(object):
    """Class keeping parameters for the em algorithm.

    """

    def __init__(self,
                 flat_prior=0.01,
                 inclusion_app=True,
                 inclusion_list=None,
                 max_iteration=15,
                 gaussian_norm=True,
                 t1c_extra_hyper=True):
        """Constructor.

        Parameters
        ----------
        flat prior: np.float [0.01]
        inclusion_app: boolean [True]
            Whether to include inter-modality constraints.
        inclusion_list: None or list of list of str or int [None]
            If None, default inclusion constraints will be used,
            If list, specifies a list of pairs of modalities,
            in which the tumor voxels of the first modality are
            constraint to be a strict subset of the tumor voxels
            of the second modality.
            Modality should be specified by string ['t1', 't1c', 't2',
            'flair'] or by integer [0,1,2,3].
        max_iteration: int [15]
            how many iterations the em has to do.
        guassian_norm: boolean [True]
            whether the gaussian evaluations have to be normalised.
        t1c_extra_hyper: boolean [True]
            when t1c hyperintensities are expected much more intense
            than the healthy tissue intensities.

        """
        self.flat_prior = flat_prior
        self.inclusion_app = inclusion_app
        self.set_inclusion_list(inclusion_list)
        self.max_iteration = max_iteration
        self.gaussian_norm = gaussian_norm
        self.t1c_extra_hyper = t1c_extra_hyper

    def set_inclusion_list(self, inclusion_list):
        """Set the inclusion list as a list of pairs of integers or a list
        of pairs of strings, where the integers refer to the indices of
        the data_files supplied in ReaderEM(), or the strings directly
        to modalities.

        parameters
        ----------
         - inclusion_list: None or
                          list of pairs of strings or integers [None]
            If None, default inclusion constraints will be used,
            If list, specifies a list of pairs of modalities,
            in which the tumor voxels of the first modality are
            constraint to be a strict subset of the tumor voxels
            of the second modality.
            Modality should be specified by string ['t1', 't1c', 't2',
            'flair'] or by integer [0,1,2,3].
        """

        if inclusion_list != [] and \
                inclusion_list is not None:
            if True in [len(pair) != 2 for pair in inclusion_list]:
                raise AttributeError('inclusion_list in Param(),' + \
                                     'should contain PAIRS of ' + \
                                     'strings or integers')
            for pair in inclusion_list:
                for elem in pair:
                    if not (isinstance(elem, int) or isinstance(elem, str)):
                        err = 'inclusion_list can only ' + \
                            'contain pairs of STRINGs or INTEGERs'
                        raise AttributeError(err)
                    if isinstance(elem, str):
                        if elem not in mr.MR_MODALITIES:
                            err = 'inclusion_list: string ' + \
                                'not recognized in modalities'
                            raise AttributeError(err)
                    if isinstance(elem, int):
                        if elem >= len(mr.MR_MODALITIES) or elem < 0:
                            err = 'inclusion_list: int ' + \
                                'shoud be in [0,1,2,3]' + \
                                '(t1, t1c, t2 and flair)'
                            raise AttributeError(err)
        self.inclusion_list = inclusion_list


########################################################################


class TumorLocation(object):
    """Class to keep information about outer tumor coordinatees.

    """

    def __init__(self,
                 tumor_lower_co=None,
                 tumor_upper_co=None):

        self.tumor_lower_co = tumor_lower_co
        self.tumor_upper_co = tumor_upper_co


########################################################################


class Initializer(object):
    """Class integrating information from ems_read_input.py and
    ems_write_output.py to generate a valid initializer for
    ems_tumor_original.py.

    """

    def __init__(self, read_instance, write_instance):
        """Connects reading and writing instance to create
        valid em initializer.

        Parameters
        -----------
        read_instance: ems_read_input.ReaderEM
        write_instance: ems_write_output.WriteOuput

        """

        if not isinstance(write_instance, wo.WriterEM):
            raise ValueError('tumor_location_instance ' + \
                             'is not a TumorLocation() instance')
        self.w = write_instance
        if not isinstance(read_instance, ri.ReaderEM):
            raise ValueError('tumor_location_instance ' + \
                             'is not a TumorLocation() instance')
        self.r = read_instance

        self.w.connect_to_read(self.r)

        if self.w.abort:
            print 'Invalid write instance: ' + \
                  'objects to write exist already'
            self.invalid = True
        else:
            self.invalid = False
            self.p = self.w.r.p

    def get_metadata_path(self):
        """Get the path where metadata should be written."""

        return self.w.metadata_path

    def is_valid_for_start(self):
        """Ask if object is valid to use for ems_tumor_original.IteraterEM."""

        return not self.invalid

    def initialize(self):
        """Furter initialize read instance."""

        if not self.r.is_read:
            self.r.read()
            self.r.initialise_tumor_related()
