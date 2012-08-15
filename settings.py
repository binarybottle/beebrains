#!/usr/bin/python

"""
Settings for preprocess_beebrains.py and process_beebrains.py

Requirements for preprocess_beebrains.py:
* Python libraries:  nibabel, numpy, scipy
* ANTS registration software for motion correction
* ImageMagick -- only if converting images and creating montage movies
Requirements for process_beebrains.py:
* Python libraries:  numpy, nipy, nibabel

(c) 2012  Arno Klein, under Apache License Version 2.0
          arno@binarybottle.com  .  www.binarybottle.com
"""

import os

#-------------------------------------
# Settings for preprocess_beebrains.py
#-------------------------------------

ANTS = os.environ['ANTSPATH']  # for registration (motion correction)
IMAGEMAGICK = '/usr/local/bin/'  # only if converting images and creating montage movies

# Settings
xdim = 130  # x dimension for each image
ydim = 172  # y dimension for each image
images_per_run = 232
wavelengths = ['340', '380']
smooth_sigma = 1
ref_image = 1
ext = '.nii.gz'

# Steps to run:
convert_images = 1  # convert .pst image slice stack to 2D nifti files
correct_motion = 1  # apply registration to correct for motion
smooth_images = 0  # smooth the resulting motion-corrected images
stack_slices = 0
plot_design_matrix = 0
plot_histogram = 0
plot_contrast = 0

# Table parameters (start from 0 index):
start_row = 1
behavior_column = 1
amplitude_column = 3
wavelength_column = 4
image_file_column = 5
start1_column = 6
stop1_column = 7
start2_column = 8
stop2_column = 9

# Save image types:
save_affine = 0  #  Save affine transform -- NOTE: Using average of affine transforms instead
save_nonlinear = 0  #  Save nonlinear transform -- NOTE: Using average of nonlinear transforms instead
save_affine_avg = 1  # Use average of affine transforms
save_nonlinear_avg = 1  # Use average of nonlinear transforms
save_montages = 0  # Save for visualizing results
save_movie = 0  # Very memory intensive -- better to use another program such as GraphicConverter
