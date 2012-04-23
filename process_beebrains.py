#!/usr/bin/python

"""
Processing steps:

(1) 

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

# http://neuroimaging.scipy.org/doc/manual/html/api/generated/nipy.neurospin.utils.design_matrix.html
# class nipy.neurospin.utils.design_matrix.DesignMatrix(frametimes=None, paradigm=None, hrf_model='Canonical', drift_model='Cosine', hfcut=128, drift_order=1, fir_delays=[, 0], fir_duration=1.0, cond_ids=None, add_regs=None, add_reg_names=None)

paradigm = BlockParadigm(con_id=['awake','odor','odor','odor',...], onset=[0, 5], duration=[10, 2], amplitude=[1, 0.2, 0.3])

dmtx = make_dmtx(np.linspace(0,19,20), paradigm, hrf_model='FIR', drift_model='Blank', hfcut=np.inf)

dmtx.matrix


"""
Examples of design matrices specification and and computation
(event-related design, FIR design, etc)

Author : Bertrand Thirion: 2009-2010
https://github.com/nipy/nipy/blob/master/examples/labs/demo_dmtx.py
"""
import numpy as np
import pylab as mp
from nipy.modalities.fmri.design_matrix import make_dmtx
from nipy.modalities.fmri.experimental_paradigm import \
    EventRelatedParadigm, BlockParadigm

# frame times
tr = 1.0
nscans = 128
frametimes = np.linspace(0, (nscans - 1) * tr, nscans)

# experimental paradigm
conditions = ['c0', 'c0', 'c0', 'c1', 'c1', 'c1', 'c3', 'c3', 'c3']
onsets = [30, 70, 100, 10, 30, 90, 30, 40, 60]
hrf_model = 'Canonical'
motion = np.cumsum(np.random.randn(128, 6), 0)
add_reg_names = ['tx', 'ty', 'tz', 'rx', 'ry', 'rz']

# block design matrix
duration = 7 * np.ones(9)
paradigm = BlockParadigm(con_id=conditions, onset=onsets,
                             duration=duration)

X2 = make_dmtx(frametimes, paradigm, drift_model='Polynomial',
                         drift_order=3)
# plot the results
fig = mp.figure(figsize=(10, 6))
ax = mp.subplot(1, 3, 1)
X2.show(ax=ax)
ax.set_title('Block design matrix', fontsize=12)
ax = mp.subplot(1, 3, 3)


"""
This is an example where
1. An sequence of fMRI volumes are simulated
2. A design matrix describing all the effects related to the data is computed
3. A GLM is applied to all voxels
4. A contrast image is created

Author : Bertrand Thirion, 2010
https://github.com/nipy/nipy/blob/master/examples/labs/example_glm.py
"""
print __doc__

import numpy as np
import os.path as op
import tempfile

from nibabel import save, Nifti1Image
import nipy.modalities.fmri.design_matrix as dm
from nipy.labs.utils.simul_multisubject_fmri_dataset import \
     surrogate_4d_dataset
import nipy.labs.glm as GLM
from nipy.modalities.fmri.experimental_paradigm import EventRelatedParadigm

#######################################
# Simulation parameters
#######################################

# volume mask
shape = (20, 20, 20)
affine = np.eye(4)

# timing
n_scans = 128
tr = 2.4

# paradigm
frametimes = np.linspace(0, (n_scans - 1) * tr, n_scans)
conditions = np.arange(20) % 2
onsets = np.linspace(5, (n_scans - 1) * tr - 10, 20) # in seconds
hrf_model = 'Canonical'
motion = np.cumsum(np.random.randn(n_scans, 6), 0)
add_reg_names = ['tx', 'ty', 'tz', 'rx', 'ry', 'rz']

# write directory
swd = tempfile.mkdtemp()

########################################
# Design matrix
########################################

paradigm = EventRelatedParadigm(conditions, onsets)
X, names = dm.dmtx_light(frametimes, paradigm, drift_model='Cosine', hfcut=128,
                         hrf_model=hrf_model, add_regs=motion,
                         add_reg_names=add_reg_names)


#######################################
# Get the FMRI data
#######################################

fmri_data = surrogate_4d_dataset(shape=shape, n_scans=n_scans)[0]

# if you want to save it as an image
data_file = op.join(swd, 'fmri_data.nii')
save(fmri_data, data_file)

########################################
# Perform a GLM analysis
########################################

# GLM fit
Y = fmri_data.get_data()
model = "ar1"
method = "kalman"
glm = GLM.glm()
glm.fit(Y.T, X, method=method, model=model)

# specify the contrast [1 -1 0 ..]
contrast = np.zeros(X.shape[1])
contrast[0] = 1
contrast[1] = - 1
my_contrast = glm.contrast(contrast)

# compute the constrast image related to it
zvals = my_contrast.zscore()
contrast_image = Nifti1Image(np.reshape(zvals, shape), affine)

# if you want to save the contrast as an image
contrast_path = op.join(swd, 'zmap.nii')
save(contrast_image, contrast_path)


print 'wrote the some of the results as images in directory %s' % swd

h, c = np.histogram(zvals, 100)
import pylab
pylab.figure()
pylab.plot(c[: - 1], h)
pylab.title(' Histogram of the z-values')
pylab.show()