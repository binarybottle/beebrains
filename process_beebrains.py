#!/usr/bin/python

"""
Processing steps:

(1) Compute a design matrix
(2) Apply a GLM to all voxels
(3) Create a contrast image

Outputs: 

Command: python <this file name> <data directory> <output directory>

Example: python process_beebrains.py data output

Requirements:  
* Python libraries:  nipy

After:
https://github.com/nipy/nipy/blob/master/examples/labs/demo_dmtx.py
https://github.com/nipy/nipy/blob/master/examples/labs/example_glm.py

(c) 2012  Arno Klein, under Apache License Version 2.0
          arno@binarybottle.com  .  www.binarybottle.com
"""

from nipy.modalities.fmri.design_matrix import make_dmtx
from nipy.modalities.fmri.experimental_paradigm import BlockParadigm
import numpy as np
import pylab as mp
import os.path as op
import nibabel as nib
import nipy.labs.glm as GLM

ANTS = '/software/ANTS_1.9/bin/'  # for registration (motion correction)

xdim = 130  # x dimension for each image
ydim = 172  # y dimension for each image
frames_per_run = 232  # number of images captured for each run

########################################
# Construct a design matrix
########################################

frametimes = np.linspace(0,2319,2320)
contrasts = ['awake', 'odor','odor','odor','odor','odor','odor','odor','odor']
onsets = [1160, 233,465,697,929,1393,1625,1857,2089]
durations = [1160, 232,232,232,232,232,232,232,232]
amplitudes = [1, 0.6,0.4,0.3,0.2,0.6,0.4,0.3,0.2]
paradigm = BlockParadigm(con_id=contrasts, onset=onsets, duration=durations, amplitude=amplitudes)

dmtx = make_dmtx(frametimes, paradigm, hrf_model='FIR', drift_model='Blank', hfcut=np.inf)

dmtx.matrix

# plot the results
mp.figure(figsize=(10, 6))
dmtx.show()
mp.title('Block design matrix')
mp.show()

########################################
# Perform a GLM analysis
########################################

# Load data
inputvolumes = 'output/bk120309e.lst_smoothed/image*.nii.gz'
outputvolume = 'output/bk120309e.lst_smoothed.nii.gz'

cmd = [ANTS+'StackSlices', outputvolume, str(xdim), str(ydim), str(frames_per_run), inputvolumes]
print(' '.join(cmd)); os.system(' '.join(cmd))                

img = nib.load(outputvolume)
Y = img.get_data()

# GLM fit
method = "kalman"
model = "ar1"
glm = GLM.glm()
glm.fit(Y.T, X, method=method, model=model)

# specify the contrast [1 -1 0 ..]
contrast = np.zeros(X.shape[1])
contrast[0] = 1
contrast[1] = - 1
my_contrast = glm.contrast(contrast)

# compute the contrast image related to it
zvals = my_contrast.zscore()
contrast_image = nib.Nifti1Image(np.reshape(zvals, shape), affine)

# if you want to save the contrast as an image
contrast_path = op.join(swd, 'zmap.nii')
nib.save(contrast_image, contrast_path)

print 'wrote some of the results as images in directory %s' % swd

h, c = np.histogram(zvals, 100)
mp.figure()
mp.plot(c[: - 1], h)
mp.title('Histogram of the z-values')
mp.show()

