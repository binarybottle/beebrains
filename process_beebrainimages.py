#!/usr/bin/python

"""
Preprocessing steps:

(1) Open a bee's table.
(2) Convert each .pst image file listed in each row of the table
    to nifti (neuroimaging file) format.
(3) Compute an affine and a nonlinear transform for each image
    that registers that image to the middle image in its wavelength
    set (ANTS software from UPenn).
(4) For each pair of images (corresponding to wavelength set 1 and 2),
    apply their affine and nonlinear transforms.
    Also compose the average of their affine and nonlinear transforms
    and apply them to both images.
(5) Divide each motion-corrected image corresponding to one
    wavelength by the motion-corrected image of a second wavelength.
(6) Smooth each ratio image (from 5) with a Gaussian kernel.

Processing steps:

(1) Construct a design matrix
(2) Apply a GLM to all voxels
(3) Create a contrast image

Outputs: Nifti files for each table (for each bee).

Command:
python <this file> <table file> <image directory> <output directory>

Example:
python preprocess_beebrains.py data/Bee1_lr120313l.txt data/Bee1_lr120313l.pst output

Requirements:
* Python libraries:  nibabel, numpy, scipy, nipy
* ANTS registration software for motion correction
* ImageMagick -- only if creating montages/movies

fMRI-based analysis with much-appreciated help from Satrajit S. Ghosh
and after Bertrand Thirion's examples:
https://github.com/nipy/nipy/blob/master/examples/labs/demo_dmtx.py
https://github.com/nipy/nipy/blob/master/examples/labs/example_glm.py

(c) 2012  Arno Klein, under Apache License Version 2.0
          arno@binarybottle.com  .  www.binarybottle.com
"""

#-----------------------------------------------------------------------------
# Import Python libraries
#-----------------------------------------------------------------------------
import os, sys
import csv
import nibabel as nib
import numpy as np
from scipy.ndimage.filters import gaussian_filter

import pylab as mp
from nipy.modalities.fmri.design_matrix import make_dmtx
from nipy.modalities.fmri.experimental_paradigm import BlockParadigm
import nipy.labs.glm as GLM

from settings import *

#-----------------------------------------------------------------------------
# Command-line arguments
#-----------------------------------------------------------------------------
args = sys.argv[1:]
if len(args)<3:
    print("\n\t Please provide the names of two directories: \
                one containing .lst table files, another to save output.")
    print("\t Example: python " + sys.argv[0] + \
          " data/Bee1_lr1203131.txt data/Bee1_lr1203131.pst output")
    sys.exit()
else:
    table_file = str(args[0])
    images_dir = str(args[1])
    output_path = str(args[2])
    try:
        if not os.path.exists(output_path):
            os.mkdir(output_path)
    except IOError:
        print("Cannot make " + output_path + " directory.")

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------
def convert2png(stem):
    """Convert nifti file to jpeg and png"""
    cmd = [ANTS+'ConvertToJpg', stem + ext, stem + '.jpg']
    print(' '.join(cmd)); os.system(' '.join(cmd))
    cmd = [IMAGEMAGICK+'convert', stem + '.jpg', stem + '.png']
    print(' '.join(cmd)); os.system(' '.join(cmd))

def smooth2d(image_in, image_out, smooth_sigma=smooth_sigma):
    """Smooth each image with a Gaussian filter"""
    image_nib = nib.load(image_in)
    image = image_nib.get_data()
    image_smooth = gaussian_filter(image, sigma=smooth_sigma, order=0)
    image_smooth_nib = nib.Nifti1Image(image_smooth, np.eye(4))
    image_smooth_nib.to_filename(image_out)

def try_mkdir(dir_name):
    try:
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
    except IOError:
        print("Cannot make " + dir_name + " directory.")

#=============================================================================
# Read table and construct lists of values and names for output folders
#=============================================================================
try:
    csv_reader = csv.reader(open(table_file, 'rU'), dialect=csv.excel_tab)
except IOError:
    print("Cannot open " + table_file + ".")

# Loop through rows (starting from start_row) and build lists
conditions = []
amplitudes = []
onsets = []
durations = []
for irow, row in enumerate(csv_reader):
    if irow >= start_row:
        if row[behavior_column] == "wake":
            conditions.append(1)
        elif row[behavior_column] == "sleep":
            conditions.append(0)
        amplitudes.append(np.float(row[amplitude_column]))
        onsets.append(int(row[start1_column]))
        onsets.append(int(row[start2_column]))
        durations.append(int(row[stop1_column]) - int(row[start1_column]))
        durations.append(int(row[stop2_column]) - int(row[start2_column]))

# Output directory names
in_stem = os.path.splitext(os.path.basename(table_file))[0]
out_path = os.path.join(output_path, in_stem)
out_path_images = os.path.join(out_path, 'images')
out_path_transforms = os.path.join(out_path, 'transforms')
out_path_transformed = os.path.join(out_path, 'transformed')
out_path_smoothed = os.path.join(out_path, 'smoothed' + str(smooth_sigma))
out_path_montages = os.path.join(out_path, 'montages')
image_stack_file = os.path.join(out_path, 'preprocessed_slicestack' + ext)
try_mkdir(out_path)

affine = np.eye(4)
n_images = n_runs * images_per_run

#=============================================================================
# Preprocess (coregister, smooth) images
#=============================================================================
if preprocess_images:

    w1 = wavelengths[0]
    w2 = wavelengths[1]

    #-------------------------------------------------------------------------
    # Convert images to nifti files
    #-------------------------------------------------------------------------
    if convert_images:

        # Make output directory
        try_mkdir(out_path_images)

        # Load table
        try:
            csv_reader = csv.reader(open(table_file, 'rU'), dialect=csv.excel_tab)
        except IOError:
            print("Cannot open " + table_file + ".")

        # Loop through rows (starting from start_row)
        n_rows = 0
        n_runs = 0
        for irow, row in enumerate(csv_reader):
            if irow >= start_row:
                n_rows += 1
                if np.mod(n_rows, 2):
                    n_runs += 1

                # Load .pst file
                file = os.path.join(images_dir, row[image_file_column])
                wavelength = row[wavelength_column]
                print('Loading ' + file + ' and converting images...')
                raw = np.fromfile(file, dtype='<i2')

                # Loop through images
                for iframe in range(images_per_run):
                    image_vector = raw[iframe * xdim * ydim : (iframe + 1) * xdim * ydim]
                    image_matrix = np.reshape(image_vector, (xdim, ydim))

                    # Save each image as a nifti file
                    converted_file = 'run' + str(n_runs) + '_' +\
                                     'image' + str(iframe + 1) + '_' +\
                                     wavelength + ext
                    img_ratio_nib = nib.Nifti1Image(image_matrix, np.eye(4))
                    img_ratio_nib.to_filename(os.path.join(out_path_images, converted_file))

    #-------------------------------------------------------------------------
    # Apply motion correction (ANTS) to one of the images
    #-------------------------------------------------------------------------
    # Warp: ANTS 2 -m CC[target.nii.gz, source.nii.gz, 1, 2] -o transform.nii.gz
    #              -r Gauss[2,0] -t SyN[0.5] -i 30x99x11 --use-Histogram-Matching
    #              --number-of-affine-iterations 10000x10000x10000x10000x10000
    # Reslice: WarpImageMultiTransform 2 source.nii.gz source2target.nii.gz
    #              -R target.nii transformWarp.nii.gz transformAffine.txt
    #-------------------------------------------------------------------------
    if correct_motion:
        print('Correcting motion...')
        try_mkdir(out_path_transforms)
        try_mkdir(out_path_transformed)
        for irun in range(n_runs):
            for iframe in range(images_per_run):

                stem = 'run' + str(irun + 1) + '_' +\
                       'image' + str(iframe + 1)
                ref_stem = 'run' + str(irun + 1) + '_' +\
                           'image' + str(ref_image)
                image_stem = os.path.join(out_path_images, stem)
                image_ref_stem = os.path.join(out_path_images, ref_stem)
                transform_stem = os.path.join(out_path_transforms, stem)
                transformed_stem = os.path.join(out_path_transformed, stem)

                # Only run if output doesn't already exist
                if not os.path.exists(transformed_stem + '_AvgWarped_ratio' + ext):
                    for ilambda in range(len(wavelengths)):
                        w = '_' + wavelengths[ilambda]
                        image_ref = image_ref_stem + w + ext
                        image_lambda = image_stem + w + ext
                        xfm = transform_stem + w + '_' + ext
                        affine_iterations = '10000x10000x10000x10000x10000'
                        nonlin_iterations = '30x99x11'
                        cmd = [ANTS+'ANTS 2 -m CC[', image_ref + ',' ,
                               image_lambda + ',1,2] -o', xfm,
                               '-r Gauss[2,0] -t SyN[0.5] -i',
                               nonlin_iterations,
                               '--number-of-affine-iterations',
                               affine_iterations,
                               '--use-Histogram-Matching']
                    print(' '.join(cmd)); os.system(' '.join(cmd))

                    # Compute the average of the two files' affine transforms,
                    # then apply it to each of two lambda files
                    if save_affine_avg or save_nonlinear_avg:
                        cmd = [ANTS+'AverageAffineTransform 2',
                               transform_stem + '_AvgAffine.txt',
                               transform_stem + '_' + w1 + '_Affine.txt',
                               transform_stem + '_' + w2 + '_Affine.txt']
                        print(' '.join(cmd)); os.system(' '.join(cmd))
                    if save_affine_avg:
                        for ilambda in range(len(wavelengths)):
                            w = '_' + wavelengths[ilambda]
                            cmd = [ANTS+'WarpImageMultiTransform 2',
                                   image_stem + w + ext,
                                   transformed_stem + w + '_AvgAffined' + ext, '-R',
                                   image_ref_stem + w + ext,
                                   transform_stem + '_AvgAffine.txt']
                            print(' '.join(cmd)); os.system(' '.join(cmd))

                        # Divide the two average-affine motion-corrected images
                        print('Dividing average-affine motion-corrected ' +\
                              'images for each of two wavelengths...')
                        cmd = [ANTS+'ImageMath 2',
                               transformed_stem + '_AvgAffined_ratio' + ext,'/',
                               transformed_stem + '_' + w1 + '_AvgAffined' + ext,
                               transformed_stem + '_' + w2 + '_AvgAffined' + ext]
                        print(' '.join(cmd)); os.system(' '.join(cmd))

                    # Compute the average of the two files' nonlinear transforms,
                    # then apply them to the ratio of the two lambda files
                    if save_nonlinear_avg:
                        cmd = [ANTS+'AverageImages 2',
                               transform_stem + '_AvgWarp' + ext, '0',
                               transform_stem + '_*_Warp' + ext]
                        print(' '.join(cmd)); os.system(' '.join(cmd))

                        for ilambda in range(len(wavelengths)):
                            w = '_' + wavelengths[ilambda]
                            cmd = [ANTS+'WarpImageMultiTransform 2',
                                   image_stem + w + ext,
                                   transformed_stem + w + '_AvgWarped' + ext, '-R',
                                   image_ref_stem + w + ext,
                                   transform_stem + '_AvgWarp' + ext,
                                   transform_stem + '_AvgAffine.txt']
                            print(' '.join(cmd)); os.system(' '.join(cmd))

                        # Divide the two average-warp motion-corrected images
                        print('Dividing average-warp motion-corrected ' +\
                              'images for each of two wavelengths...')
                        cmd = [ANTS + 'ImageMath 2',
                               transformed_stem + '_AvgWarped_ratio' + ext,'/',
                               transformed_stem + '_' + w1 + '_AvgWarped' + ext,
                               transformed_stem + '_' + w2 + '_AvgWarped' + ext]
                        print(' '.join(cmd)); os.system(' '.join(cmd))

    #-------------------------------------------------------------------------
    # Smooth each motion-corrected image file.
    #-------------------------------------------------------------------------
    if smooth_images:
        try_mkdir(out_path_smoothed)
        for irun in range(n_runs):
            for iframe in range(images_per_run):
                print('Smoothing image ' + str(iframe + 1) +\
                      ' from run ' + str(irun + 1))
                stem = 'run' + str(irun + 1) + '_' +\
                       'image' + str(iframe + 1)
                transformed_stem = os.path.join(out_path_transformed, stem)
                smoothed_stem = os.path.join(out_path_smoothed, stem)
                if save_nonlinear_avg:
                    smooth2d(transformed_stem + '_AvgWarped_ratio' + ext,
                        smoothed_stem + '_AvgWarped_ratio_smooth' + ext,
                        smooth_sigma = smooth_sigma)

    #-------------------------------------------------------------------------
    # Convert each nifti image file to jpg and to png
    # and save a montage of preprocessed images
    #-------------------------------------------------------------------------
    if save_montages:
        try_mkdir(out_path_montages)
        for irun in range(n_runs):
            for iframe in range(images_per_run):
                stem = 'run' + str(irun + 1) + '_' +\
                       'image' + str(iframe + 1)
                print('Convert to jpg, png images; create montage for image ' + \
                      str(iframe + 1) + ' from run ' + str(irun + 1))
                image_stem = os.path.join(out_path_images, stem)
                transformed_stem = os.path.join(out_path_transformed, stem)
                smoothed_stem = os.path.join(out_path_smoothed, stem)
                out_montage = os.path.join(out_path_montages, stem + '.png')

                convert2png(image_stem + '_' + w1)
                convert2png(image_stem + '_' + w2)
                if save_affine_avg:
                    convert2png(transformed_stem + '_AvgAffined_ratio')
                if save_nonlinear_avg:
                    convert2png(transformed_stem + '_AvgWarped_ratio')
                    convert2png(smoothed_stem + '_AvgWarped_ratio_smooth')

                if save_affine_avg:
                    cmd = [IMAGEMAGICK+'montage',
                           image_stem + '_' + w1 + '.png',
                           image_stem + '_' + w2 + '.png',
                           transformed_stem + '_AvgAffined_ratio.png',
                           transformed_stem + '_AvgWarped_ratio.png',
                           smoothed_stem + '_AvgWarped_ratio_smooth.png',
                           '-label ' + str(iframe+1),
                           '-geometry +1+0 -tile 5x -background black',
                           out_montage]
                else:
                    cmd = [IMAGEMAGICK+'montage',
                           image_stem + '_' + w1 + '.png',
                           image_stem + '_' + w2 + '.png',
                           transformed_stem + '_AvgWarped_ratio.png',
                           smoothed_stem + '_AvgWarped_ratio_smooth.png',
                           '-label ' + str(iframe+1),
                           '-geometry +1+0 -tile 4x -background black',
                           out_montage]
                print(' '.join(cmd)); os.system(' '.join(cmd))

    #-------------------------------------------------------------------------
    # Stack preprocessed images for analysis
    #-------------------------------------------------------------------------
    if stack_slices:
        print('Stack preprocessed images...')
        shape = (xdim, ydim, 1, n_images)
        image_stack = np.zeros(shape)
        for irun in range(n_runs):
            for iframe in range(images_per_run):
                stem = 'run' + str(irun + 1) + '_' + \
                       'image' + str(iframe + 1)
                smoothed = os.path.join(out_path_smoothed,
                                        stem + '_AvgWarped_ratio_smooth' + ext)
                img = nib.load(smoothed)
                slice = img.get_data()
                image_stack[:,:,0,i-1] = slice
        img = nib.Nifti1Image(image_stack, affine)
        nib.save(img, image_stack_file)

#=============================================================================
# Conduct a general linear model analysis on the preprocessed images
# (Requires image_stack and the following paradigm lists from above:
#  conditions, onsets, durations, amplitudes)
#=============================================================================
if analyze_images:
    if not preprocess_images and not stack_slices:
        img = nib.load(image_stack_file)
        image_stack = img.get_data()

    #-------------------------------------------------------------------------
    # Construct a design matrix
    #-------------------------------------------------------------------------
    frametimes = np.linspace(0, n_images-1, n_images)

    paradigm = BlockParadigm(con_id=conditions, onset=onsets,
                             duration=durations, amplitude=amplitudes)

    dmtx = make_dmtx(frametimes, paradigm, hrf_model='FIR',
                     drift_model='Blank', hfcut=np.inf)
    design_matrix = dmtx.matrix

    # Plot the design matrix
    if plot_design_matrix:
        fig1 = mp.figure(figsize=(10, 6))
        dmtx.show()
        mp.title('Block design matrix')
        fig1_file = os.path.join(out_path, 'design_matrix.png')
        mp.savefig(fig1_file)

    #-------------------------------------------------------------------------
    # Apply a general linear model to all pixels
    #-------------------------------------------------------------------------
    method = "kalman"
    model = "ar1"
    glm = GLM.glm()
    glm.fit(image_stack.T, design_matrix, method=method, model=model)

    #-------------------------------------------------------------------------
    # Create a contrast image
    #-------------------------------------------------------------------------
    # Specify the contrast [1 -1 0 ..]
    contrast = np.zeros(design_matrix.shape[1])
    contrast[0] = 1
    #contrast[1] = 0
    #contrast[2] = -1
    my_contrast = glm.contrast(contrast)

    # Compute the contrast image
    affine = np.eye(4)
    zvals = my_contrast.zscore()
    contrast_image = nib.Nifti1Image(zvals.T, affine)

    # Save the contrast as an image
    contrast_file = os.path.join(out_path, 'zmap' + ext)
    nib.save(contrast_image, contrast_file)

    # Plot histogram
    if plot_histogram:
        h, c = np.histogram(zvals, 100)
        fig2 = mp.figure()
        mp.plot(c[: - 1], h)
        mp.title('Histogram of the z-values')
        fig2_file = os.path.join(out_path, 'histogram.png')
        mp.savefig(fig2_file)

    # Plot contrast image
    if plot_contrast:
        fig3 = mp.figure()
        mp.matshow(zvals)
        mp.title('Contrast image')
        fig3_file = os.path.join(out_path, 'contrast.png')
        mp.savefig(fig3_file)
