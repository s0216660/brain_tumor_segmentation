# -*- coding: UTF-8 -*-
"""Module allowing to write data and metadata resulting from the em
algorithm.

"""
__version__ = '0.2'
__author__ = 'Esther Alberts'

import os, shutil
import pickle

from utils import own_itk as oitk

import ems_initializer as ei

########################################################################


class WriterEM(object):
    """Class to deal with (over)writing data and results from the algorithm.

    An instance must be created to create an ems_initializer.Initializer()
    instance, which is needed to initialize the
    ems_tumor_original.IteraterEM(), which allows to run the em.

    Attributes
    ----------
    overwrite : boolean
        If True, results are to be overwritten if they already exist,
        if False, abort is set to True.
    abort : boolean
        If True, the em should not be allowed to run.
    param_string : str
        Suffix to append to paths where tissue and tumor segmentations
        will be written.
    save_dir : str
        Path to directory were tissue and tumor segmentations are to be
        written.
    tissue_paths : list of str
        Paths where tissue segmentations will be saved.
    tumor_paths : list of str
        Paths where tumor segmentations will be saved.
    metadata_path : str
        Path were metadata is to be written.
    metadata : dict
        Dictionary containing metadata information.
    r : ems_read_input.ReaderEM
        Reader containing information about input data and run parameters.

    """
    def __init__(self,
                 overwrite=False,
                 param_string='',
                 save_dir=None):
        """Constructor.

        Parameters
        ----------
        overwrite : boolean, optional
            If True, results are to be overwritten if they already exist,
            if False, abort is set to True.
        param_string : str, optional
            Suffix to append to paths where tissue and tumor segmentations
            will be written.
        save_dir : str, optional
            Path to directory were tissue and tumor segmentations are to be
            written.

        """

        # parameters to write files
        self.overwrite = overwrite
        self.param_string = param_string
        self.save_dir = save_dir
        self.r = None

    ####################################################################

    def connect_to_read(self, read_instance):
        """Integrate information from reading instance.

        Parameters
        ----------
        read_instance : ems_read_input.ReadInput
            the read_instance allows to have access to the input data.

        """

        print 'Initializing write instance with read instance'
        self.r = read_instance

        # Further refine WriterEM instance
        self._set_directories()
        self._create_param_file()
        self._set_result_paths()
        self._check_existence()

    ####################################################################

    def _set_directories(self):
        """Set the directory where results will be saved."""
        # Create directory to store results

        if self.save_dir is None:
            self.save_dir = os.path.join(\
                os.path.dirname(self.r.data_files[0]),
                'ems_results')

        if not os.path.exists(self.save_dir):
            os.mkdir(self.save_dir)

    def _create_param_file(self):
        """Write a parameter file and keep a metadata.pickle instance"""

        self.metadata_path = os.path.join(self.save_dir,
                                          'metadata.dat')
        self.tmp_metadata_path = os.path.join(self.save_dir,
                                              'metadata_tmp.dat')
        fid = open(self.tmp_metadata_path, 'w')
        self.metadata = {}

        fid.write('// \n')
        fid.write('// ****** Input Data ****** \n')
        fid.write('atlas_files (gm, wm, csf): %s \n' \
                  % str(self.r.atlas_files))
        fid.write('mask_files (gm, wm, csf): %s \n' \
                  % str(self.r.mask_files))
        fid.write('hyper_data_files : %s \n' \
                  % str(self.r.hyper_data_files))
        fid.write('hypo_data_files : %s \n' \
                  % str(self.r.hypo_data_files))
        fid.write('data_files : %s \n' % str(self.r.data_files))
        fid.write('modalities : %s \n' % str(self.r.modalities))
        input_data = {}
        input_data['atlas_files'] = self.r.atlas_files
        input_data['mask_files'] = self.r.mask_files
        input_data['hyper_data_files'] = self.r.hyper_data_files
        input_data['hypo_data_files'] = self.r.hypo_data_files
        input_data['data_files'] = self.r.data_files
        input_data['modalities'] = self.r.modalities
        self.metadata['input_data'] = input_data

        fid.write('// ****** Tumor Box ****** \n')
        fid.write('tumor_lower_co : %s \n' % str(self.r.tlc))
        fid.write('tumor_upper_co : %s \n' % str(self.r.tuc))
        tumor_box = {}
        tumor_box['tumor_lower_co'] = self.r.tlc
        tumor_box['tumor_upper_co'] = self.r.tuc
        self.metadata['tumor_box'] = tumor_box

        fid.write('// ****** Parameters ****** \n')
        fid.write('inclusion : %s \n' \
                  % str(self.r.inclusion))
        fid.write('inclusion list [t1, t1c, t2, flair] : %s \n' \
                  % str(self.r.inclusion_list))
        fid.write('max_iteration : %s \n' \
                  % str(self.r.p.max_iteration))
        fid.write('gauss_per_class : %s \n' \
                  % str(self.r.gauss_per_class))
        fid.write('flat_prior : %s \n' \
                  % str(self.r.p.flat_prior))
        fid.write('t1c_extra_hyper : %s \n' \
                  % str(self.r.p.t1c_extra_hyper))
        param = {}
        param['inclusion'] = self.r.inclusion
        param['inclusion_list'] = self.r.inclusion_list
        param['max_iteration'] = self.r.p.max_iteration
        param['gauss_per_class'] = self.r.gauss_per_class
        param['flat_prior'] = self.r.p.flat_prior
        if self.r.p.t1c_extra_hyper:
            param['t1c_extra_hyper'] = self.r.p.t1c_extra_hyper
        self.metadata['param'] = param

        fid.write('// ****** Extra parameters ****** \n')
        fid.write('gaussian_norm : %s \n' \
                  % str(self.r.p.gaussian_norm))
        extra_param = {}
        extra_param['gaussian_norm'] = self.r.p.gaussian_norm
        self.metadata['extra_param'] = extra_param

        fid.close()

    def _set_result_paths(self):
        """Set the paths where tissue and tumor segmentations should
        be written.

        """

        tissue_paths = []
        for i in range(3):
            tissue_paths.append(os.path.join(self.save_dir,
                ei.TISSUES[i] + self.param_string + '.nii'))

        tumor_paths = []
        for i in range(len(self.r.data_files)):
            mod = self.r.modalities[i]
            tumor_paths.append(os.path.join(self.save_dir,
                mod + '_tumor' + self.param_string + '.nii'))

        self.tissue_paths = tissue_paths
        self.tumor_paths = tumor_paths


    def _check_existence(self):
        """Check if the paths to write results already exist."""

        self.abort = False
        tumor_paths_exist = [os.path.exists(path) \
                             for path in self.tumor_paths]
        tissue_paths_exist = [os.path.exists(path) \
                             for path in self.tissue_paths]
        if not (False in tumor_paths_exist) and \
                not (False in tissue_paths_exist) and \
                    os.path.exists(self.metadata_path):
            print 'This data seems to have already been ' + \
                  'processed (possibly with different parameters)'

            if not self.overwrite:
                self.abort = True
            else:
                print 'Pre-calculated files are NOW being removed !!!'
                for the_file in os.listdir(self.save_dir):
                    file_path = os.path.join(self.save_dir, the_file)
                    if file_path != self.tmp_metadata_path:
                        try:
                            if os.path.isfile(file_path):
                                os.remove(file_path)
                            elif os.path.isdir(file_path):
                                shutil.rmtree(file_path)
                        except Exception, err:
                            print err

        if not self.abort:
            pickle_path = os.path.join(self.save_dir,
                                       'metadata.pickle')
            pickle.dump(self.metadata,
                        open(pickle_path, 'w'))

            os.rename(self.tmp_metadata_path,
                      self.metadata_path)

    ####################################################################

    def write_tissue_maps(self, tissue_maps):
        """Write the tissue segmentations."""

        for class_nr in range(self.r.nr_classes):

            segment_tissue_class = \
                self.r.to_matrix(tissue_maps[class_nr, :])
            image = oitk.make_itk_image(
                    segment_tissue_class, self.r.data_itk_image)

            if class_nr == ei.GM:
                oitk.write_itk_image(image, self.tissue_paths[0])
            if class_nr == ei.WM:
                oitk.write_itk_image(image, self.tissue_paths[1])
            if class_nr == ei.CSF:
                oitk.write_itk_image(image, self.tissue_paths[2])

        print 'Tissue segmentation written to ' + self.tissue_paths[0]
        print 'Tissue segmentation written to ' + self.tissue_paths[1]
        print 'Tissue segmentation written to ' + self.tissue_paths[2]

    def write_tumor_maps(self, tumor_maps, tumor_atlas):
        """Write the tumor segmentations."""

        # Write tumor segmentation.
        tumor_atlas_full = \
            self.r.to_matrix(tumor_atlas[ei.TUMOR, :])
        tumor_image = oitk.make_itk_image(tumor_atlas_full,
                                          self.r.data_itk_image)
        result_file = os.path.join(self.save_dir,
                                   'tumor_atlas.nii')
        oitk.write_itk_image(tumor_image, result_file)

        # Write channel-specific tumor segmentations.
        for channel in range(self.r.nr_channels):
            tumor_channel = \
                self.r.to_matrix(tumor_maps[channel, :])
            tumor_image = oitk.make_itk_image(
                tumor_channel, self.r.data_itk_image)
            oitk.write_itk_image(tumor_image,
                                 self.tumor_paths[channel])
            print 'Tumor segmentation written to '\
                  + self.tumor_paths[channel]

    ####################################################################
