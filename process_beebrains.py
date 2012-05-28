#!/usr/bin/python

"""
Processing steps:

(1) Construct a design matrix
(2) Apply a GLM to all voxels
(3) Create a contrast image

Outputs: 

Command: python <this file name> <subject string> <data directory>

Example: python process_beebrains.py bk120309e.lst output

Requirements:  
* Python libraries:  numpy, nipy, nibabel

With help from Satrajit S. Ghosh and after Bertrand Thirion's examples:
https://github.com/nipy/nipy/blob/master/examples/labs/demo_dmtx.py
https://github.com/nipy/nipy/blob/master/examples/labs/example_glm.py

(c) 2012  Arno Klein, under Apache License Version 2.0
          arno@binarybottle.com  .  www.binarybottle.com
"""

import os, sys
import numpy as np
import pylab as mp
import nibabel as nib
from nipy.modalities.fmri.design_matrix import make_dmtx
from nipy.modalities.fmri.experimental_paradigm import BlockParadigm
import nipy.labs.glm as GLM

########################################
# Settings
########################################

xdim = 130  # x dimension for each image
ydim = 172  # y dimension for each image
frames_per_run = 232  # number of images captured for each run
num_runs = 10

# table name:
# subject = 'bk120309e.lst'  
# conditions -- the bee woke up or an odor concentration was introduced:
conditions = ['awake', 'odor','odor','odor','odor','odor','odor','odor','odor']  
# when these conditions began, by frame number in concatenated sequence:
onsets = [1160, 233,465,697,929,1393,1625,1857,2089]
# durations of each condition:
durations = [1160, 232,232,232,232,232,232,232,232]
# the amplitude for each condition (awake = 1; concentration < 1):
amplitudes = [1, 0.6,0.4,0.3,0.2,0.6,0.4,0.3,0.2]

stack_slices = 0
plot_design_matrix = 0
plot_histogram = 0
plot_contrast = 0

smooth_sigma = 1

# Command-line arguments
args = sys.argv[1:]
if len(args)<2:
    print("\n\t Please provide a subject string and the name of the data/output directory.")
    print("\t Example: python " + sys.argv[0] + " bk120309e.lst output")
    sys.exit()
else:
    subject = str(args[0]) # 'bk120309e.lst'
    out_path = str(args[1])

preprocessed_path = os.path.join(out_path, subject+'_smoothed' + str(smooth_sigma) +'/')
preprocessed_volume = os.path.join(out_path, subject+'_smoothed' + str(smooth_sigma) + '.nii.gz')

########################################
# Construct a design matrix
########################################

nframes = frames_per_run * num_runs
frametimes = np.linspace(0,nframes-1,nframes)

paradigm = BlockParadigm(con_id=conditions, onset=onsets, duration=durations, amplitude=amplitudes)

dmtx = make_dmtx(frametimes, paradigm, hrf_model='FIR', drift_model='Blank', hfcut=np.inf)
X = dmtx.matrix

# Plot the design matrix
if plot_design_matrix:
    fig1 = mp.figure(figsize=(10, 6))
    dmtx.show()
    mp.title('Block design matrix for ' + subject)
    fig1_file = os.path.join(out_path, subject + '_design_matrix.png')
    mp.savefig(fig1_file)
    #mp.show()

########################################
# Apply a GLM to all voxels
########################################

shape = (xdim,ydim,1,nframes)
affine = np.eye(4)

# Load data
if stack_slices:
    Y = np.zeros(shape)
    for i in range(1,nframes+1):       
        img = nib.load(os.path.join(preprocessed_path, 
                       'image'+str(i)+'_AvgWarped_ratio_smooth.nii.gz'))
        slice = img.get_data()
        Y[:,:,0,i-1] = slice
    img = nib.Nifti1Image(Y, affine)
    nib.save(img, preprocessed_volume)
else:
    img = nib.load(preprocessed_volume)
    Y = img.get_data()

# GLM fit
method = "kalman"
model = "ar1"
glm = GLM.glm()
glm.fit(Y.T, X, method=method, model=model)

########################################
# Create a contrast image
########################################

# Specify the contrast [1 -1 0 ..]
contrast = np.zeros(X.shape[1])
contrast[0] = 1
#contrast[1] = 0
#contrast[2] = -1
my_contrast = glm.contrast(contrast)

# Compute the contrast image
zvals = my_contrast.zscore()
contrast_image = nib.Nifti1Image(zvals.T, affine)

# Save the contrast as an image
contrast_file = os.path.join(out_path, subject+'_zmap_smooth' + str(smooth_sigma) + '.nii.gz')
nib.save(contrast_image, contrast_file)

# Plot histogram
if plot_histogram:
    h, c = np.histogram(zvals, 100)
    fig2 = mp.figure()
    mp.plot(c[: - 1], h)
    mp.title('Histogram of the z-values for ' + subject)
    fig2_file = os.path.join(out_path, subject + '_histogram_sigma' + str(smooth_sigma) + '.png')
    mp.savefig(fig2_file)
    #mp.show()

# Plot contrast image
if plot_contrast:
    fig3 = mp.figure()
    mp.matshow(zvals)
    mp.title('Contrast image for ' + subject)
    fig3_file = os.path.join(out_path, subject + '_contrast_sigma' + str(smooth_sigma) + '.png')
    mp.savefig(fig3_file)
    #mp.show()

