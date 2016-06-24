# -*- coding: UTF-8 -*-
"""Module responsible for loading images contained in image files in numpy
arrays and setting variables needed to run the em on.

"""
__version__ = '0.2'
__author__ = 'Esther Alberts'

import numpy as np

from utils import own_itk as oitk
from utils import mr_variables as mr
from utils import normalisation as norm

import ems_initializer as ei

########################################################################


class InputData(object):
    """Class to read input data and initialise variables based on
    run parameters.

    An instance must be created to create an ems_initializer.Initializer()
    instance, which is needed to initialize the
    ems_tumor_original.IteraterEM(), which allows to run the em.

    Attributes
    ----------
    p : ems_initializer.Param
        parameter instance.
    tuc : ndarray, int
        upper, right, back coordinates of the tumor box.
    tlc : ndarray, int
        lower, left, front coordinates of the tumor box.
    nr_channels : int
        number of modalities.
    data_files : 1D array, str
        path to each input modality.
    modalities : 1D array, str
        name of the modality corresponding to each of the paths in data_files.
    hyper : 1D array, boolean
        True if the corresponding path in data_files shows a hyperintense
        tumor, False if hypointense tumor.
    data_itk_image : itk image
        Proto itk image to be used when results have to be converted from
        ndarrays to itk images.
    dim : tuple of int
        shape of the input image ndarrays.
    data_spaces : tuple of float
        pixel spacing of the input images.
    data_type : data_type
        data_type of the input images.
    full_data : 2D array
        contains the flattened input images along the first axis.
    inclusion : boolean
        True if inclusion constraints are to be applied across the modalities.
    inclusion_list : list of list of int
        List of pairs where each pair of integers.

    """

    def __init__(self,
                 hyper_data_files,
                 hypo_data_files,
                 param_instance,
                 tumor_location):
        """Constructor.

        Parameters:
        -----------
        hyper_data_files : str, list of str or dict of (str, str)
            Paths to images in which tumor appears hyperintense, if
            dict, keys are in ['t1c', 'flair', 't2'].
        hypo_data_files : str, list of str or dict of (str, str)
            Paths to images in which tumor appears hypointense, if
            dict, keys are in ['t1'].
        param_instance : ems_initializer.Param()
            Parameter instance.
        tumor_location : ems_initializer.TumorLocation()
            Tumor location instance.

        .. note:: make sure that all the images referenced in
            hyper_data_files and hypo_data_files are in the same
            reference space and have the same pixel resolution and
            image dimensions.

        """

        self._set_data_files(hyper_data_files, hypo_data_files)

        if not isinstance(param_instance, ei.Param):
            raise ValueError('param_instance is not a Param() instance')
        else:
            self.p = param_instance
        self.set_inclusion(self.p.inclusion_app, self.p.inclusion_list)

        if tumor_location is None:
            tumor_location = ei.TumorLocation()
        if not isinstance(tumor_location, ei.TumorLocation):
            err = '`tumor_location_instance` is no `TumorLocation()` instance'
            raise ValueError(err)
        else:
            if tumor_location.tumor_lower_co is None:
                self.tlc = np.zeros(3, dtype='int')
            else:
                self.tlc = tumor_location.tumor_lower_co
            if tumor_location.tumor_upper_co is None:
                self.tuc = self.dim - np.ones(3, dtype='int')
            else:
                self.tuc = tumor_location.tumor_upper_co


    ####################################################################

    def _set_data_files(self, hyper_data_files, hypo_data_files):
        """Define paths to input brain tumor images.

        Set self.data_itk_image, self.dim, self.data_spaces,
        self.data_type and self.full_data.

        """
        # Concatenate into _data_files and set related instance variables
        self._concatenate_data_files(hyper_data_files, hypo_data_files)

        # 1. Check data dimensionalities:
        #     spacing, dim and data_type should be the same
        # 2. Read data in a numpy array
        dim = 0
        data_spaces = 0
        data_type = 0
        for channel in range(self.nr_channels):
            print 'Reading image at ' + self.data_files[channel] + '...'
            image, tmp_dim, tmp_spaces = \
                    oitk.get_itk_data(self.data_files[channel], True)
            if channel == 0:
                data = np.zeros(shape=(self.nr_channels,
                                       np.prod(tmp_dim)),
                                dtype=image.dtype)
            else:
                if tmp_dim != dim or tmp_spaces != data_spaces \
                        or image.dtype != data_type:
                    err = 'Please make sure data dimensions ' + \
                          'and pixel spacing are consistent ' + \
                          'along all channels'
                    raise ValueError(err)
            flat_image = np.ravel(image)
            data[channel, :] = flat_image
            dim = tmp_dim
            data_spaces = tmp_spaces
            data_type = image.dtype

        self.data_itk_image = oitk.get_itk_image(self.data_files[0])
        self.dim = dim
        self.data_spaces = data_spaces
        self.data_type = data_type

        self.full_data = data

    def _concatenate_data_files(self, hyper_data_files, hypo_data_files):
        """Define paths to input brain tumor images.

        Set self.modalities, self.data_files, self.nr_channels and
        self.hyper.

        """
        if hyper_data_files is None and hypo_data_files is None:
            raise AttributeError('None objects as input')
        if hyper_data_files is None:
            hyper_data_files = []
        if hypo_data_files is None:
            hypo_data_files = []
        if isinstance(hyper_data_files, str):
            hyper_data_files = [hyper_data_files]
        if isinstance(hypo_data_files, str):
            hypo_data_files = [hypo_data_files]
        if isinstance(hyper_data_files, dict):
            these_keys = hyper_data_files.keys()
            these_values = hyper_data_files
            hyper_mod = []
            hyper_data_files = []
            for key in mr.HYPER_MODALITIES:
                if key in these_keys:
                    hyper_mod.append(key)
                    hyper_data_files.append(these_values[key])
            if len(these_keys) > len(hyper_mod):
                raise ValueError('Modality key in hyper_data_files' + \
                                 ' not recognized or not hyper intense')
        else:
            hyper_mod = mr.get_modality_strings(hyper_data_files)
        if isinstance(hypo_data_files, dict):
            these_keys = hypo_data_files.keys()
            these_values = hypo_data_files
            hypo_mod = []
            hypo_data_files = []
            for key in mr.HYPO_MODALITIES:
                if key in these_keys:
                    hypo_mod.append(key)
                    hypo_data_files.append(these_values[key])
            if len(these_keys) > len(hypo_mod):
                raise ValueError('Modality key in hypo_data_files' + \
                                 ' not recognized or not hypo intense')
        else:
            hypo_mod = mr.get_modality_strings(hypo_data_files)

        self.hyper_data_files = hyper_data_files
        self.hypo_data_files = hypo_data_files

        data_files = np.concatenate((hyper_data_files,
                                     hypo_data_files),
                                    axis=0)
        modalities = np.concatenate((hyper_mod,
                                     hypo_mod),
                                    axis=0)
        oitk.check_path_existence(data_files)
        nr_channels = len(data_files)
        hyper = np.zeros(shape=(nr_channels,))
        hyper[0:len(hyper_data_files)] = 1

        self.data_files = data_files
        self.modalities = modalities
        self.nr_channels = nr_channels
        self.hyper = hyper

    def is_t1c_channel(self, channel):
        """Return whether an index belongs to the t1c channel."""
        return self.modalities[channel] == mr.T1C

    ####################################################################

    def set_inclusion(self, inclusion_app, inclusion_list=None):
        """Set wheather inclusion constraints have to be activated,
        and if so, specify which inclusion constraints should hold
        (indexed by mr_variables.MR_MODALITIES)."""

        inclusion = True
        if np.any(self.modalities == '?'):
            inclusion = False
            print('Inclusion will NOT hold')
        if inclusion_app == False:
            inclusion = False
            print('Inclusion will NOT hold')

        if inclusion:
            if inclusion_list is None:
                self.inclusion_list = [[0, 3], [1, 3], [2, 3]]
            else:
                self.inclusion_list = inclusion_list
        else:
            self.inclusion_list = []

        self.inclusion = inclusion

    def get_inclusion(self):
        """Get the list of inclusion pairs."""

        if hasattr(self, 'inclusion'):
            return self.inclusion
        else:
            raise RuntimeError('Please call set_inclusion first')


########################################################################


class ReaderEM(InputData):
    """Class to read input data and initialise variables based on
    run parameters.

    Attributes
    ----------
    p : ems_initializer.Param
        parameter instance.
    tuc : ndarray, int
        upper, right, back coordinates of the tumor box.
    tlc : ndarray, int
        lower, left, front coordinates of the tumor box.
    nr_channels : int
        number of modalities.
    data_files : 1D array, str
        path to each input modality.
    modalities : 1D array, str
        name of the modality corresponding to each of the paths in data_files.
    hyper : 1D array, boolean
        True if the corresponding path in data_files shows a hyperintense
        tumor, False if hypointense tumor.
    data_itk_image : itk image
        Proto itk image to be used when results have to be converted from
        ndarrays to itk images.
    dim : tuple of int
        shape of the input image ndarrays.
    data_spaces : tuple of float
        pixel spacing of the input images.
    data_type : data_type
        data_type of the input images.
    full_data : 2D array
        contains the flattened input images along the first axis.
    inclusion : boolean
        True if inclusion constraints are to be applied across the modalities.
    inclusion_list : list of list of int
        Specifies a list of pairs of modalities, in which the tumor
        voxels of the first modality are constraint to be a strict subset
        of the tumor voxels of the second modality.
        Integers are referring to the modalities in modalities.
    mask_files: None, str or list of str
        strings indicate filenames containing the binary mask(s)
        of one or more of the supplied images.
    atlas_files: list of str or dict (str, str)
        Paths of GM, WM and CSF registered tissue atlases
        If dict, keys should be ['gm', 'wm','csf'] and values strings,
        else the order of grey matter, white matter and cerebrospinal
        fluid *respectively* should be respected.
    gauss_per_class: list of int
        number of gaussians used to model GM, WM and CSF respectivily.
    nr_of_classes : int
        Number of tissues to be modelled, is constant and set to three
        (WM, GM and CSF).
    nr_gaussians : int
        Number of guassians present to model all tissue classes.
    lkp : list of int
        Specifies for every gaussian which tissue class it is used to model.
    playing : 1D array
        Specifies for every voxel whether it is a brain voxel.
    nr_brain_voxels : int
        Number of brain voxels.
    atlas : 2D array
        Registered tissue atlases (nr_casses, nr_brain_voxels).
    data : 2D array
        Input images (nr_channels, nr_brain_voxels).
    atlas_prior : 2D array
        Registered tissue atlases per gaussian
        (nb_gaussians, nr_brain_voxels).
    gauss_maps_init : 2D array
        Set to atlas_prior.
    tumor_mask : 1D array
        Specifies for every brain voxel whether it is in the tumor box.
    tumor_atlas_init : 2D array
        Specifies for every brain voxel how likeli it is to be tumor
        (first axis) or how likeli it is not to be tumor (second axis).
    tumor_maps_init : 2D array
        Initialisation of tumor segmentation for brain voxels
        (nr_channels, nr_brain_voxels).
    nr_tissue_conf : int
        number of tissue-tumor configuration accross the channels.
    tumor_transitions : 2D array
        tumor-non-tumor configurations accross the channels for every
        tissue-tumor configuration.
    classes : 1D array
        Specifies for every tissue-tumor configuration which tissue is
        being considered.
    gaussians : 1D array
        Specifies for every tissue-tumor configuration which gaussian is
        being considered.

    """

    def __init__(self, hyper_data_files,
                       hypo_data_files,
                       mask_files,
                       atlas_files,
                       param_instance,
                       tumor_location=None,
                       gauss_per_class=None):
        """Read all data and allocate in variables.

        Parameters
        ----------
        hyper_data_files : str, list of str or dict of (str, str)
            Paths to images in which tumor appears hyperintense, if
            dict, keys are in ['t1c', 'flair', 't2'].
        hypo_data_files : str, list of str or dict of (str, str)
            Paths to images in which tumor appears hypointense, if
            dict, keys are in ['t1'].
        mask_files: None, str or list of str
            strings indicate filenames containing the binary mask(s)
            of one or more of the supplied images.
        atlas_files: list of str or dict (str, str)
            Paths of GM, WM and CSF registered tissue atlases
            If dict, keys should be ['gm', 'wm','csf'] and values strings,
            else the order of grey matter, white matter and cerebrospinal
            fluid *respectively* should be respected.
        param_instance : ems_initializer.Param()
            Parameter instance.
        tumor_location : ems_initializer.TumorLocation()
            Tumor location instance.
        gauss_per_class: list of int
            number of gaussians used to model GM, WM and CSF respectivily.

        .. NOTE:: All files should contain preprocessed images, in the
            sense that, all images (referenced in hyper_data_files,
            hypo_data_files, mask_files and atlas_files) should have the
            same resolution, the same isotropic pixel dimensions and
            they should be in the same reference space.
            Accepted file extensions are: *TIFF, JPEG, PNG, BMP,
            DICOM, GIPL, Bio-Rad, LSM, Nifti, Analyze,
            SDT/SPR (Stimulate), Nrrd or VTK images*.

        """

        InputData.__init__(self, hyper_data_files,
                                 hypo_data_files,
                                 param_instance,
                                 tumor_location)

        self.mask_files = mask_files
        self.atlas_files = atlas_files

        if gauss_per_class is None:
            gauss_per_class = [1, 1, 1]
        self.gauss_per_class = gauss_per_class

        self.is_read = False

    ####################################################################

    def read(self):
        """Function reading in all data and setting most of the
        important instance variables.

        """

        print 'Reading the data, loading images, ' + \
                'initializing numpy ndarrays'

        raw_mask = self._set_mask_files()
        raw_atlas = self._set_atlas_files()

        self._set_mask(raw_atlas, raw_mask)
        self._set_atlas(raw_atlas)
        self._set_data()

        self._set_tissue_gaussians()
        self._initialise_gauss_maps()

        self.is_read = True

    ####################################################################

    def _set_mask_files(self):
        """Function reading the mask files.

        """
        # Read brain mask(s).
        if self.mask_files != None:
            oitk.check_path_existence(self.mask_files)
            playing = np.ones(shape=(np.prod(self.dim)), dtype='?')
            if isinstance(self.mask_files, str):
                self.mask_files = [self.mask_files]
            else:
                print('You supplied ' + str(len(self.mask_files)) + ' masks')
                print('The intersection will be taken as the true mask')
            for file_nb in range(len(self.mask_files)):
                mask, dim_mask = oitk.get_itk_data(\
                                        self.mask_files[file_nb])[0:-1]
                if dim_mask != self.dim:
                    raise ValueError('Please make sure mask dimensions ' + \
                        'correspond to data dimensions')
                playing = np.logical_and(playing, np.ravel(mask))
        else:
            playing = None

        return playing

    def _set_atlas_files(self):
        """Function reading the atlas files."""

        if not isinstance(self.atlas_files, list):
            if isinstance(self.atlas_files, dict):
                if set(self.atlas_files.keys()) != \
                            set(['gm', 'wm', 'csf']):
                    raise AttributeError('keys in dict atlas_files ' + \
                                         'should be "gm", "wm" and "csf"')
                atlas_files = [self.atlas_files['gm'],
                               self.atlas_files['wm'],
                               self.atlas_files['csf']]
                self.atlas_files = atlas_files
            else:
                raise AttributeError('atlas_files should be list or dict')
        self.nr_classes = 3

        # Read atlas.
        if len(atlas_files) != self.nr_classes:
            err = 'Atlas files should be of lenght 3: WM, GM and CSF'
            raise AttributeError(err)

        atlas = np.zeros(shape=(self.nr_classes, np.prod(self.dim)))
        for class_nr in range(self.nr_classes):
            atlas_raw, dim_atlas = \
                oitk.get_itk_data(atlas_files[class_nr])[0:-1]
            if dim_atlas != self.dim:
                raise ValueError('Please make sure atlas dimensions ' + \
                        'correspond to data dimensions')
            atlas[class_nr, :] = np.ravel(atlas_raw)
        atlas = norm.window(atlas, 0, 1)

        return atlas

    ####################################################################

    def _set_mask(self, atlas, playing=None):
        """Set the mask depending on the atlas.

        Set self.playing, self.nr_brain_voxels.

        """
        sum_atlas = np.sum(atlas, axis=0)

        if playing is None:
            # Define brain voxels in case mask is not given
            playing = np.ones(shape=(np.prod(self.dim),), dtype='?')
            playing[sum_atlas > 0.5] = True
            playing[sum_atlas <= 0.5] = False
            self.playing = playing
        else:
            # Mask is given, now throw out brain voxels with nonzero
            # atlas probabilities.
            self.playing = np.logical_and(sum_atlas > 0, playing)

        self.nr_brain_voxels = np.count_nonzero(self.playing)

    def _set_atlas(self, atlas):
        """Normalise, flatten and mask the atlas.

        Set self.atlas.

        """
        sum_atlas = np.sum(atlas, axis=0)

        self.atlas = np.zeros(shape=(self.nr_classes,
                                     self.nr_brain_voxels))
        for class_nr in range(self.nr_classes):
            atlas[class_nr, self.playing] = \
                    atlas[class_nr, self.playing] / sum_atlas[self.playing]
            self.atlas[class_nr, :] = atlas[class_nr, self.playing]

    def _set_data(self):
        """Flatten and mask the data.

        Set self.data.

        """
        self.data = np.zeros(shape=(self.nr_channels,
                                    self.nr_brain_voxels),
                             dtype=self.data_type)
        for channel in range(self.nr_channels):
            channel_data = self.full_data[channel, :]
            self.data[channel, :] = channel_data[self.playing]

    ####################################################################

    def _set_tissue_gaussians(self):
        """Specify gaussians for each tissue class

        Set self.lkp, self.nr_gaussians.

        """

        lkp = []  # Class to which a certain gaussian belongs
        for i in range(len(self.gauss_per_class)):
            lkp.extend(np.ones(self.gauss_per_class[i], dtype='int') * (i))
        self.lkp = np.array(lkp)
        self.nr_gaussians = len(lkp)

    def _initialise_gauss_maps(self):
        """Initialise guassian-specific tissue classifications and set
        the atlas prior.

        Set self.gauss_maps_init, self.atlas_prior.

        """
        # Initialise healthy tissue classification.
        tissue_maps = np.zeros(shape=(self.nr_classes,
                                      self.nr_brain_voxels))
        atlas_prior = np.zeros(shape=(self.nr_gaussians,
                                      self.nr_brain_voxels))
        for gauss in range(self.nr_gaussians):
            class_nr = self.lkp[gauss]
            tissue_maps[class_nr, :] = self.atlas[class_nr, :]
            atlas_prior[gauss, :] = self.atlas[class_nr, :] / \
                                    self.gauss_per_class[class_nr]

        self.atlas_prior = atlas_prior
        self.gauss_tissue_maps_init = atlas_prior

    ####################################################################
    # Functions to be called after read() has been executed
    ####################################################################

    def initialise_tumor_related(self):
        """Initialise tumor classification.

        Initialise all tumor-related variables, such as tumor mask,
        tumor atlas, and tumor configurations across the channels.

        """
        if not self.is_read:
            err = 'Function cannot be called before this instance read ' + \
                  'out the data first (call self.read())'
            raise RuntimeError(err)

        self._set_tumor_box()
        self._set_tumor_atlas()
        self._set_tumor_configurations()

    ####################################################################

    def _set_tumor_box(self):
        """Create tumor mask.

        Set self.tumor_mask.

        """
        # Create tumor mask
        print 'original image dimensions: ' + str(self.dim)

        dim_box = (self.tuc - self.tlc) + np.ones(3, dtype='int')
        print 'tumor box dimensions: ' + str(dim_box)
        tumor_mask = np.zeros(self.dim)
        tumor_mask[self.tlc[0]:self.tuc[0] + 1,
                   self.tlc[1]:self.tuc[1] + 1,
                   self.tlc[2]:self.tuc[2] + 1] = 1
        tumor_mask = np.ravel(tumor_mask)[self.playing]
        self.tumor_mask = tumor_mask

    def _set_tumor_atlas(self):
        """Initialise tumor atlas.

        Set self.tumor_atlas_init.

        """
        # Initialise tumor atlas
        tumor_atlas = np.zeros(shape=(2, self.nr_brain_voxels))
        tumor_atlas[ei.TUMOR, :] = self.tumor_mask * self.p.flat_prior
        tumor_atlas[ei.NO_TUMOR, :] = 1 - tumor_atlas[ei.TUMOR, :]
        self.tumor_atlas_init = tumor_atlas

        self._initialise_tumor_maps()


    def _initialise_tumor_maps(self):
        """Initialise tumor classification.

        Set self.tumor_maps_init.

        """
        # Initialise tumor classification.
        tumor_maps = np.zeros(shape=(self.nr_channels,
                                     self.nr_brain_voxels))
        for channel in range(self.nr_channels):
            tumor_maps[channel, :] = self.tumor_atlas_init[ei.TUMOR, :]
        self.tumor_maps_init = tumor_maps

    def _set_tumor_configurations(self):
        """Create tumor configurations.

        Set self.nr_tissue_conf, self.tumor_transitions
        self.classes, self.gaussians.

        """
        # Specifiy tumor classes: tumor is allowed only in WM and GM
        tumor_classes = np.zeros((self.nr_classes,), dtype='?')
        tumor_classes[ei.WM] = 1
        tumor_classes[ei.GM] = 1
        no_tumor_classes = 1 - tumor_classes

        # Specify tumor presence configurations along the channels.
        nrTumorConf = 2 ** self.nr_channels
        tumor_transition_matrix = np.zeros(shape=(nrTumorConf,
                                                  self.nr_channels),
                                           dtype='?')
        for conf in range(nrTumorConf):
            tumor_transition_matrix[conf, :] = \
                [bool(int(i)) for i in \
                    list(format(conf, '0' + str(self.nr_channels) + 'b'))]

        # throw out unrealistic tumor transitions
        if self.get_inclusion():
            tumor_transition_matrix = \
                self.__enforce_inclusion(tumor_transition_matrix)
        nr_tumor_trans = tumor_transition_matrix.shape[0]

        # Calculate matrix with all tumor configurations
        # within all tissue classes.
        nr_tissue_conf = np.dot(self.gauss_per_class, no_tumor_classes) + \
                         (np.dot(self.gauss_per_class, tumor_classes) * \
                          nr_tumor_trans)
        classes = np.zeros(shape=(nr_tissue_conf))
        gaussians = np.zeros(shape=(nr_tissue_conf))
        tumor_transitions = np.zeros(shape=(nr_tissue_conf,
                                            self.nr_channels),
                                     dtype='?')
        ind = 0
        for gauss in range(self.nr_gaussians):
            class_nr = self.lkp[gauss]
            if tumor_classes[class_nr] == ei.NO_TUMOR:
                classes[ind] = class_nr
                gaussians[ind] = gauss
                tumor_transitions[ind, :] = 0
                ind = ind + 1
            if tumor_classes[class_nr] == ei.TUMOR:
                classes[ind:ind + nr_tumor_trans] = class_nr
                gaussians[ind:ind + nr_tumor_trans] = gauss
                tumor_transitions[ind:ind + nr_tumor_trans, :] = \
                    tumor_transition_matrix
                ind = ind + nr_tumor_trans

        self.nr_tissue_conf = nr_tissue_conf
        self.tumor_transitions = tumor_transitions
        self.classes = classes
        self.gaussians = gaussians

    def __enforce_inclusion(self, matrix):
        """Remove certain tumor configurations from the configuration
        matrix. E.g., tumor cannot be present in t1c if it is not
        present in flair.

        """

        all_modalities = mr.MR_MODALITIES

        indices = []
        for i in range(len(all_modalities)):
            this_mod = all_modalities[i]
            this_ind = np.where(np.asarray(self.modalities) == this_mod)[0]
            if len(this_ind) == 1:
                this_ind = this_ind[0]
            else:
                this_ind = None
            indices.append(this_ind)

        nr_combinations = 0
        for i in range(matrix.shape[0]):
            # get a row from the input matrix
            combination = np.zeros([1, matrix.shape[1]], dtype='?')
            combination[0, :] = matrix[i, :]

            # check if the row should be considered or not
            log = 0
            for pair in self.inclusion_list:
                if not isinstance(pair[0], int):
                    small_ind = np.where(np.asarray(all_modalities) \
                                         == pair[0])[0]
                else:
                    small_ind = pair[0]
                if not isinstance(pair[1], int):
                    big_ind = np.where(np.asarray(all_modalities) \
                                       == pair[1])[0]
                else:
                    big_ind = pair[1]
                if not None in [indices[small_ind], indices[big_ind]]:
                    log += (combination[0][indices[small_ind]] \
                                > combination[0][indices[big_ind]])

            # add the row if it should be considered
            if log == 0:
                nr_combinations = nr_combinations + 1
                if nr_combinations == 1:
                    restricted_matrix = combination
                else:
                    restricted_matrix = np.concatenate((restricted_matrix,
                                                        combination))

        return restricted_matrix

    ####################################################################

    def to_matrix(self, flat_array):
        """Rescale flat array(s) back into 3D volumes.

        Rescale flat array, or the array of flat arrays, back into 3D
        volumes, corresponding to the shape of the input images. Each
        flat array should should either be of length brain voxels,
        or of the length of all voxels present in the 3d volume.

        """

        def is_playing(arr):
            """Return whether the array contains only the brain voxels
            or all voxels of the entire 3d volume """
            if arr.size == self.nr_brain_voxels:
                return True
            elif arr.size == np.prod(dim):
                return False
            else:
                err = 'this array does not have an interpreteble size'
                raise RuntimeError(err)

        dim = self.dim

        if len(flat_array.shape) == 1:
            if is_playing(flat_array):
                flat_array_full = np.zeros((np.prod(dim)),)
                flat_array_full[self.playing] = flat_array
            else:
                flat_array_full = flat_array

            matrix = flat_array_full.reshape(dim)

        elif len(flat_array.shape) == 2:
            flat_arrays = flat_array
            shape = (flat_array.shape[0], dim[0], dim[1], dim[2])
            matrix = np.zeros(shape=shape)

            for i in range(flat_array.shape[0]):
                flat_array = flat_arrays[i]
                matrix[i, :, :, :] = self.__to_matrix(flat_array)

        return matrix

    ####################################################################
