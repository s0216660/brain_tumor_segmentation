# -*- coding: UTF-8 -*-
"""Module containing functions enabling to read, make and
write ITK images.

"""
__version__ = '0.2'
__author__ = 'Esther Alberts'

import os as os
import SimpleITK as itk

########################################################################

def make_itk_image(arr, proto_image=None):
    """Create an itk image given an image array.
    
    Parameters
    ----------
    arr : ndarray
        Array to create an itk image with.
    proto_image : itk image, optional
        Proto itk image to provide Origin, Spacing and Direction.
        
    Returns
    -------
    image : itk image
        The itk image containing the input array `arr`.

    """

    image = itk.GetImageFromArray(arr)
    if proto_image != None:
        image.CopyInformation(proto_image)

    return image

def write_itk_image(image, path):
    """Write an itk image to a path.
    
    Parameters
    ----------
    image : itk image
        Image to be written.
    path : str
        Path where the image should be written to.

    """

    writer = itk.ImageFileWriter()
    writer.SetFileName(path)

    if os.path.splitext(path)[1] == '.nii':
        Warning('You are converting nii, ' + \
                'be careful with type conversions')

    writer.Execute(image)

def get_itk_image(path):
    """Get an itk image given a path.
    
    Parameters
    ----------
    path : str
        Path pointing to an image file with extension among 
        *TIFF, JPEG, PNG, BMP, DICOM, GIPL, Bio-Rad, LSM, Nifti, Analyze,
        SDT/SPR (Stimulate), Nrrd or VTK images*.
        
    Returns
    -------
    image : itk image
        The itk image.

    """

    reader = itk.ImageFileReader()
    reader.SetFileName(path)

    image = reader.Execute()

    return image

def get_itk_array(path_or_image):
    """ Get an image array given a path or itk image.
    
    Parameters
    ----------
    path_or_image : str or itk image
        Path pointing to an image file with extension among 
        *TIFF, JPEG, PNG, BMP, DICOM, GIPL, Bio-Rad, LSM, Nifti, Analyze,
        SDT/SPR (Stimulate), Nrrd or VTK images* or an itk image.
    
    Returns
    -------
    arr : ndarray
        Image ndarray contained in the given path or the itk image.

    """

    if isinstance(path_or_image, str):
        image = get_itk_image(path_or_image)
    else:
        image = path_or_image

    arr = itk.GetArrayFromImage(image)

    return arr

def get_itk_data(path_or_image, verbose=False):
    """Get the image array, image size and pixel dimensions given an itk
    image or a path.

    Parameters
    ----------
    path_or_image : str or itk image
        Path pointing to an image file with extension among 
        *TIFF, JPEG, PNG, BMP, DICOM, GIPL, Bio-Rad, LSM, Nifti, Analyze,
        SDT/SPR (Stimulate), Nrrd or VTK images* or an itk image.
    verbose : boolean, optional
        If true, print image shape, spacing and data type of the image
        corresponding to `path_or_image.`
    
    Returns
    -------
    arr : ndarray
        Image array contained in the given path or the itk image.
    shape : tuple
        Shape of the image array contained in the given path or the itk 
        image.
    spacing : tuple
        Pixel spacing (resolution) of the image array contained in the
        given path or the itk image.
        
    """

    if isinstance(path_or_image, str):
        image = get_itk_image(path_or_image)
    else:
        image = path_or_image

    arr = itk.GetArrayFromImage(image)
    shape = arr.shape
    spacing = image.GetSpacing()[::-1]
    data_type = arr.dtype

    if verbose:
        print '\t image shape: ' + str(shape)
        print '\t image spacing: ' + str(spacing)
        print '\t image data type: ' + str(data_type)

    return arr, shape, spacing

def check_path_existence(list_of_files, raise_exception=True):
    """Check if all pasts present in list_of_files exist.

    Parameters
    ----------
    list_of_files : str, list of str or list of list of str
        every str should refer to a path, if one of the paths
        doesn't exist, False is returned.
    raise_exception : boolean, optional
        if True, an exception is raise if one of the paths doesn't
        exist.

    """

    # flatten the list
    if isinstance(list_of_files, str):
        list_of_files = [list_of_files]
    while isinstance(list_of_files[0], list):
        tmp = [item for sublist in list_of_files for item in sublist]
        list_of_files = tmp

    # check if paths exist
    for path in list_of_files:
        if not os.path.exists(path):
            if raise_exception:
                raise AttributeError('Path does not exist: ' + path)
            else:
                print('Path does not exist: ' + path)
                return False
    return True
