#!/usr/bin/python

"""
Processing steps:

(1) Construct a design matrix
(2) Apply a GLM to all voxels
(3) Create a contrast image

Example:
python preprocess_beebrains.py data/Bee1_lr120313l.txt output

Requirements:  
* Python libraries:  numpy, nipy, nibabel

With help from Satrajit S. Ghosh and after Bertrand Thirion's examples:
https://github.com/nipy/nipy/blob/master/examples/labs/demo_dmtx.py
https://github.com/nipy/nipy/blob/master/examples/labs/example_glm.py

(c) 2012  Arno Klein, under Apache License Version 2.0
          arno@binarybottle.com  .  www.binarybottle.com
"""

import pylab as mp
from nipy.modalities.fmri.design_matrix import make_dmtx
from nipy.modalities.fmri.experimental_paradigm import BlockParadigm
import nipy.labs.glm as GLM

from settings import *

include n_runs and preprocess,analyze_images in settings
add read_table in prep & n_images = n_runs * images_per_run

#=============================================================================
# Preprocess (coregister, smooth) images
#=============================================================================
if preprocess_images:

    #-------------------------------------------------------------------------
    # Stack preprocessed images
    #-------------------------------------------------------------------------
    image_stack_file = os.path.join(out_path, 'preprocessed_slicestack' + ext)
    affine = np.eye(4)

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
# (Requires image_stack and the following paradigm lists:
#  conditions, onsets, durations, amplitudes)
#=============================================================================
if analyze_images:
    if not stack_slices:
        img = nib.load(image_stack_file)
        image_stack = img.get_data()

    #-------------------------------------------------------------------------
    # Construct a design matrix
    #-------------------------------------------------------------------------
    frametimes = np.linspace(0, n_images-1, n_images)

    paradigm = BlockParadigm(con_id=conditions, onset=onsets,
                             duration=durations, amplitude=amplitudes)

    dmtx = make_dmtx(frametimes, paradigm, hrf_model='FIR', drift_model='Blank',
                     hfcut=np.inf)
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