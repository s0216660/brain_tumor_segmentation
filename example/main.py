#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""Module to run the example.

"""
__version__ = '0.2'
__author__ = 'Esther Alberts'

import os

from em import ems_initializer as ei
from em import ems_write_output as ew
from em import ems_read_input as ri
from em import ems_tumor_original as em

from utils import mr_variables as mr

########################################################################


def arrange_files():
    """Retrieve all the paths necessary to run the em."""

    # Specify data directory
    example_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(example_dir,
                            'HGG_brats_2013_pat0004_1')

    # Specify filenames of data channels
    brats_dir = os.path.join(data_dir,
                             'preprocessed_data')
    def in_brats_dir(filename):
        """Get full path of image_files filenames."""
        return os.path.join(brats_dir,
                            filename)
    image_files = {}
    filename = 'VSD.Brain.XX.O.MR_Flair.703.nii'
    image_files[mr.FLAIR] = in_brats_dir(filename)
    filename = 'VSD.Brain.XX.O.MR_Flair.703_mask.nii'
    brats_flair_mask = in_brats_dir(filename)
    filename = 'VSD.Brain.XX.O.MR_T1.704.nii'
    image_files[mr.T1] = in_brats_dir(filename)
    filename = 'VSD.Brain.XX.O.MR_T1c.705.nii'
    image_files[mr.T1C] = in_brats_dir(filename)
    filename = 'VSD.Brain.XX.O.MR_T2.706.nii'
    image_files[mr.T2] = in_brats_dir(filename)

    # Specify filenames of GM, WM and CSF atlas maps registered to image data
    atlas_dir = os.path.join(data_dir,
                             'registered_atlas')
    def in_atlas_dir(filename):
        """Get full path of atlas filenames."""
        return os.path.join(atlas_dir,
                            filename)
    atlas = {}
    atlas[ei.TISSUES[ei.GM]] = in_atlas_dir('grey_reg_via_mask.nii')
    atlas[ei.TISSUES[ei.WM]] = in_atlas_dir('white_reg_via_mask.nii')
    atlas[ei.TISSUES[ei.CSF]] = in_atlas_dir('csf_reg_via_mask.nii')

    # Store paths in variables
    # You can fill in your own filenames, but make sure all images
    # are in the same reference space with the same pixel dimensions
    hyper_files = {mod : image_files[mod]
                   for mod in mr.HYPER_MODALITIES \
                   if mod in image_files}
    hypo_files = {mod : image_files[mod]
                  for mod in mr.HYPO_MODALITIES \
                  if mod in image_files}
    mask_file = brats_flair_mask
    atlas_files = {tissue  : atlas[tissue]
                   for tissue in ei.TISSUES}

    # Direcory to save results
    save_dir = os.path.join(data_dir,
                            'ems_results')

    return {'save':save_dir,
            'hyper':hyper_files,
            'hypo':hypo_files,
            'mask':mask_file,
            'atlas':atlas_files}


def get_initializer(path_dict):
    """Create em reader and writer and get the em initializer."""

    param_instance = ei.Param(flat_prior=0.1,
                              inclusion_list=None)
    read_instance = ri.ReaderEM(path_dict['hyper'],
                                path_dict['hypo'],
                                path_dict['mask'],
                                path_dict['atlas'],
                                param_instance=param_instance)

    # Create writer
    write_instance = ew.WriterEM(overwrite=True,
                                 param_string='_example',
                                 save_dir=path_dict['save'])

    # Create initializer by connecting reader and writer
    init = ei.Initializer(read_instance=read_instance,
                          write_instance=write_instance)

    return init

def run():
    """Run the algorithm on the example."""

    # Get all relevant paths
    path_dict = arrange_files()

    # Set parameters and create reader
    init = get_initializer(path_dict)

    # Run tumor segmentation
    if init.is_valid_for_start():
        init.initialize()
        em_instance = em.IteraterEM(init)
        em_instance.run()

        # get paths of results
        tumor_paths = em_instance.w.tumor_paths
        tissue_paths = em_instance.w.tissue_paths
    else:  # paths already exist
        # get existing paths
        tumor_paths = init.w.tumor_paths
        tissue_paths = init.w.tissue_paths
        # print existing paths
        print 'em Already calculated for these paths'
        print 'Tumor segmentation paths ' + str(tumor_paths)
        print 'Tissue segmentation paths ' + str(tissue_paths)

########################################################################

if __name__ == '__main__':

    run()
