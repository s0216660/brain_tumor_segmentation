#Expectation-Maximisation algorithm for brain tumor segmentation
 * **Version**: 2016.2  
 * **Date**: 20/06/2016  
 * **Author**: Esther Alberts

---
#General usage notes

This code allows a probabilistic tumor segmentation based on [1] and [2]. Spatial regularisation by means of Markov random fields is not included.
		
##Install

Python 2.7 should be installed. 

In a terminal, run
```
python setup.py install
```
to complete the installation.

##Code execution

  * Only one command line script is provided, which will run an example. To do so, run  
```
python main.py
```
  * To use another dataset - in `main.py` set:   
    1. Path(s) to T1, T1c, T2 and/or FLAIR images in the `image_files` dictionary. Supply at least one path, at most all four paths. The algorithm will look for hypointense tumors in T1 and for hyperintense tumors in T1c, T2 and FLAIR.
    2. Paths to wm, gm and csf atlas files in the `atlas` dictionary. Please supply all three of the atlas files.
    3. Path to a brain mask in `brats_flair_mask`, optional, can be set to None.
    3. Path to store results in `save_dir`, optional, can be set to None.   
__Make sure that all your images, mask(s) and atlas maps (GM, WM and CSF) are in the same reference space and have the same resolution and dimension__!

##Results

Results are written to a directory the user specifies via `save_dir` in `ems_write_output.WriteOutput()`. By default, results are 
written in a subdirectory ems_results of the directory of the input MR images. This directory, ems_results/, contains:

> gm<*param_string*>.nii
> :  grey matter segmentation of the subjects brain
	
> wm<*param_string*>.nii
> :  white matter segmentation of the subjects brain

> csf<*param_string*>.nii
> :  csf matter segmentation of the subjects brain

> <*modality*>_tumor<*param_string*>.nii
> :  tumor segmentation for modality <modality>,

where `param_string` is set in `ems_write_output.WriteOutput()`.

##Contents

em/
> ems_read_input.py
> :  Module dealing with reading input data and initializing variables based on given parameters

> ems_write_output.py
> :  Module dealing with writing out results (to be overwritten or not) and storing metadata file containing run-specific parameters

> ems_initializer.py
> :  Module connecting information from ems_read_input.py and ems_write_output.py to generate a valid initializer for ems_tumor_original.py

> ems_tumor_original.py
> :  Module containing the iterative EM implementation

utils/
> distributions.py
> :  Module allowing to define and evaluate distributions

> own_itk.py
> :  Contains functions enabling to read, make and write ITK images

> mr_variables.py
> :  Module containing information specific to MR modalities

> normalisation.py
> :  Module including basic functions for data normalisation

example/
> HGG_brats_2013_pat0004_1/
> :  Dataset of a patient suffering from a high-grade tumor

> main.py
> :  Script including an example for tumor segmentation

---
#Contact

Bugs, problems or questions can be adressed at:

 * **E-mail**: esther.alberts@tum.de
 * **Institute**: Technical University Munchen, Klinikum rechts der Isar

---

[1] [Menze BH, Van Leemput K, Lashkari D, Weber M-A, Ayache N, Golland P. A Generative Model for Brain Tumor Segmentation in Multi-Modal Images. Medical image computing and computer-assisted intervention : International Conference on Medical Image Computing and Computer-Assisted Intervention (MICCAI) 2010; 13(Pt 2) : pp. 151-159.](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3050038/)  
[2] [Menze BH, Van Leemput K, Lashkari D, et al. A Generative Probabilistic Model and Discriminative Extensions for Brain Lesion Segmentation - With Application to Tumor and Stroke. IEEE Transactions on Medical Imaging; vol. 35, no. 4, pp. 933-946, April 2016.](http://dx.doi.org/10.1109/TMI.2015.2502596)
