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

#-----------------------------------------------------------------------------
# Settings
#-----------------------------------------------------------------------------
xdim = 130  # x dimension for each image
ydim = 172  # y dimension for each image
images_per_run = 232  # number of images for a given set of conditions (or bee)
onset_list = [73, 93]
duration_list = [11, 11]
amplitude_list = [0.000001, 0.0001, 0.001, 0.01]
wavelengths = ['340', '380']  # strings indicating paired images to divide
smooth_sigma = 1  # sigma of Gaussian kernel
ext = '.nii.gz'  # output file extension

#-----------------------------------------------------------------------------
# Run processing steps (1=True, 0=False)
#-----------------------------------------------------------------------------
convert_images = 0  # convert .pst image slice stack to 2D nifti files
divide_images = 0  # divide one wavelength's image volume by the other
correct_motion = 0  # apply registration to correct for motion
smooth_images = 0  # smooth the resulting motion-corrected images
run_analysis = 1
tests = [1, 0, 0, 0, 0]
plot_design_matrix = 1
plot_histogram = 0
plot_contrast = 1

#-----------------------------------------------------------------------------
# Table parameters (indices start from 0):
#-----------------------------------------------------------------------------
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
save_montages = 0
