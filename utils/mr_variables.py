# -*- coding: UTF-8 -*-
"""Module containing information specific to MR modalities.

"""
__version__ = '0.2'
__author__ = 'Esther Alberts'

import os
import numpy as np

########################################################################

T1 = 't1'
T1C = 't1c'
T2 = 't2'
FLAIR = 'flair'

MR_MODALITIES = [T1, T1C, T2, FLAIR]
HYPER_MODALITIES = [T1C, T2, FLAIR]
HYPO_MODALITIES = [T1]

########################################################################

def get_modality_path(files, modality):
    """Return the file in `files` which corresponds with `modality`."""

    for file_nb in range(len(files)):
        path = files[file_nb]
        mod = get_modality_string(os.path.basename(path))
        if mod == modality:
            return path, file_nb

    return 0, len(files)

def remove_modality(files, modality):
    """Remove the file in `files` which corresponds with `modality`."""

    removals = 0
    for file_nb in range(len(files)):
        path = files[file_nb - removals]
        mod = get_modality_string(path)
        if mod == modality:
            files.remove(path)
            removals += 1

    return files

def get_hyper_paths(files, modalities=None):
    """Get the files in `files` which correspond to modalities present
    in `HYPER_MODALITIES`. The modalities of the files can be specified
    with `modalities`.
    
    """

    return get_paths('hyper', files, modalities)

def get_hypo_paths(files, modalities=None):
    """Get the files in `files` which correspond to modalities present
    in `HYPO_MODALITIES`. The modalities of the files can be specified
    with `modalities`.
    
    """

    return get_paths('hypo', files, modalities)

def get_paths(selection_type, files, modalities=None):
    """ Get the files in `files` which correspond to modalities present
    in `HYPER_MODALITIES` (`selection_type` == 'hyper') or `HYPO_MODALITIES`
    (`selection_type` == 'hypo'). The modalities of the files can be
    specified with `modalities`.
    
    """

    if selection_type == 'hyper':
        type_modalities = HYPER_MODALITIES
    elif selection_type == 'hypo':
        type_modalities = HYPO_MODALITIES
    else:
        err = 'Unvalid selection type, only hyper or hypo'
        raise AttributeError(err)

    sel_files = []
    sel_mods = []
    modalities = get_modalities(files, modalities)
    for mod in type_modalities:
        for i in range(len(files)):
            if mod == modalities[i]:
                sel_files.append(files[i])
                sel_mods.append(modalities[i])

    return sel_files, sel_mods

def order_by_modality(files, modalities=None):
    """Order the files in `files` conform the order in `MR_MODALITIES`."""

    modalities = get_modalities(files, modalities)
    if len(files) > 4:
        raise ValueError('This function is written for T1, T1C, T2 ' + \
                ' and FLAIR, but you have given more than 4 files')

    new_files = [''] * 4

    for mod in MR_MODALITIES:
        for i in range(len(files)):
            if mod == modalities[i]:
                new_files.append(files[i])

    return new_files

def get_modalities(files, modalities, error=True):
    """Get the modalities of the files in `files`. Prior knowledge of the
    modalities can be in `modalities`, unknown modalities are specified by
    None.
    
    """

    if modalities is None:
        modalities = get_modality_strings(files)
    if None in modalities:
        new_modalities = get_modality_strings(files)
        modalities = [new_modalities[i] if modalities[i] is None \
                      else modalities[i] \
                      for i in range(len(modalities))]
    if valid_modalities(files, modalities, error=error):
        return modalities
    else:
        return None

def valid_modalities(files, modalities, error=True):
    """Return True if the `modalities` are valid for the specified
    `files`. The `modalities` should all be specified and be present
    in `MR_MODALITIES`."""

    msg = []
    if len(files) != len(modalities):
        msg += 'Modalities ' + str(modalities) + \
                    ' dont match files ' + str(files)

    elif False in [m in MR_MODALITIES for m in modalities]:
        msg += '\nUnvalid modality supplied ' + str(modalities)

    if msg:
        if error:
            raise ValueError(msg)
        else:
            print msg
            return False
    else:
        return True

def get_modality_string(description, interact=False):
    """Automatically derive what `modality` should be associated with
    the `description` string.
    
    """

    mod_presence = {mod:0 for mod in MR_MODALITIES}
    mod_strings = {mod:[] for mod in MR_MODALITIES}
    mod_strings[FLAIR] = ['FLAIR', 'Flair', 'FLAIR']
    mod_strings[T1] = ['T1', 'T1', 'mpr', 'MPR']
    mod_strings[T1C] = ['T1C', 'T1c', 't1C', 'T1C',
                        't1w', 'T1w', 't1W', 'T1W',
                        'KM', 'CA',
                        'gd', 'Gd', 'GD', 'ga', 'Ga', 'GA']
    mod_strings[T2] = ['T2', 'T2']

    for mod in MR_MODALITIES:
        if np.any([(string in description) \
                   for string in mod_strings[mod]]):
            mod_presence[mod] = 1
            modality = mod

    # if T1C is present, T1 might wrongly be considered to be present
    if mod_presence[T1C]:
        mod_presence[T1] = 0
        modality = T1C

    occurences = mod_presence.values().count(1)
    if occurences != 1:
        print 'Modality type is unclear, ' + \
                'check description: ' + str(description)
        if interact:
            msg = 'Do you want to continue with this file? (0/1)\n'
            valid_file = int(raw_input(msg))
            if valid_file:
                valid_mod = False
                while not valid_mod:
                    msg = 'Which modality is it? (T1,T1C,T2,FLAIR)\n'
                    modality = raw_input(msg)
                    valid_mod = (modality in MR_MODALITIES)
                    msg += 'Please type one of these: (T1,T1C,T2,FLAIR)\n'
            else:
                err = 'Unrecognized modality, user chose to quit'
                raise RuntimeError(err)
        else:
            err = 'Unrecognized modality for ' + description
            raise RuntimeError(err)

    return modality

def get_modality_strings(descriptions, interact=False):
    """Automatically derive what `modalities` should be associated with
    the `descriptions` strings."""

    modalities = [''] * len(descriptions)
    for i in range(len(descriptions)):
        description = descriptions[i]
        mod = get_modality_string(description, interact)
        modalities[i] = mod

    return modalities
