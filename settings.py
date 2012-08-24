#!/usr/bin/python

"""
Settings for process_beebrainimages.py

Requirements for process_beebrainimages.py:
* Python libraries:  nibabel, numpy, scipy, nipy
* ANTS registration software for motion correction
* ImageMagick -- only if creating montages/movies

(c) 2012  Arno Klein, under Apache License Version 2.0
          arno@binarybottle.com  .  www.binarybottle.com
"""

import os

#-----------------------------------------------------------------------------
# Settings
#-----------------------------------------------------------------------------
ANTS = os.environ['ANTSPATH']  # for registration (motion correction)
ANTSAFFINE = os.environ['ANTSPATH2']  # for registration (motion correction)
IMAGEMAGICK = ''  # only if converting images and creating montage movies

# Settings
xdim = 130  # x dimension for each image
ydim = 172  # y dimension for each image
images_per_run = 232  # number of images for a given set of conditions (or bee)
wavelengths = ['340', '380']  # strings indicating paired images to divide
smooth_sigma = 1  # sigma of Gaussian kernel
ref_image = 1  # image number to use as registration reference
ext = '.nii.gz'  # output file extension

#-----------------------------------------------------------------------------
# Run preprocessing steps (1=True, 0=False)
#-----------------------------------------------------------------------------
preprocess_images = 1
#-----------------------------------------------------------------------------
convert_images = 0  # convert .pst image slice stack to 2D nifti files
correct_motion = 1  # apply registration to correct for motion
compute_nonlinear = 0  # Compute a nonlinear transform (else just affine)
smooth_images = 1  # smooth the resulting motion-corrected images
stack_slices = 1  # save slice-stack of preprocessed images

#-----------------------------------------------------------------------------
# Run processing steps (1=True, 0=False)
#-----------------------------------------------------------------------------
analyze_images = 0
#-----------------------------------------------------------------------------
plot_design_matrix = 1
plot_histogram = 1
plot_contrast = 1

#-----------------------------------------------------------------------------
# Table parameters (indices start from 0):
#-----------------------------------------------------------------------------
start_row = 1
behavior_column = 1
amplitude_column = 3
wavelength_column = 4
image_file_column = 5
start1_column = 6
stop1_column = 7
start2_column = 8
stop2_column = 9

#-----------------------------------------------------------------------------
# Visualizing results
#-----------------------------------------------------------------------------
save_montages = 1
