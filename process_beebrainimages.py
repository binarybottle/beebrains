#!/usr/bin/python
"""
Preprocessing steps:

(1) Open a bee's table.
(2) Convert each .pst image file listed in each row of the table
    to nifti (neuroimaging file) format and create a slice stack
    corresponding to each of two wavelengths.
(3) Apply FSL's motion correction to each slice stack.
(4) Divide the motion-corrected image volume for one
    wavelength by the motion-corrected image for the second wavelength.
(6) Smooth each ratio image with a Gaussian kernel.

Processing steps:

(1) Construct a design matrix
(2) Apply a GLM to all voxels
(3) Create a contrast image

Contrast condition 1 vs. condition 2, holding condition 3 constant,
in our case, sleep vs. awake holding concentration of odorant constant.

Outputs: Nifti files for each table (for each bee).

Command:
python <this file> <table file> <image directory> <output directory>

Example:
python process_beebrainimages.py data/Bee1_lr120313l.txt data/Bee1_lr120313l.pst output

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
import nibabel as nb
import numpy as np
import pylab as mp
from nipy.modalities.fmri.design_matrix import make_dmtx
from nipy.modalities.fmri.experimental_paradigm import BlockParadigm
import nipy.labs.glm as GLM

from settings import *

#-----------------------------------------------------------------------------
# Command-line arguments and output file names
#-----------------------------------------------------------------------------
args = sys.argv[1:]
if len(args)<3:
    print("\n\t Please provide the names of two directories: \
                one containing .lst table files, another to save output.")
    print("\t Example: python " + sys.argv[0] + \
          " data/Bee1_lr120313l.txt data/Bee1_lr120313l.pst output")
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

# Output directory
in_stem = os.path.splitext(os.path.basename(table_file))[0]
out_path = os.path.join(output_path, in_stem)
try:
    if not os.path.exists(out_path):
        os.mkdir(out_path)
except IOError:
    print("Cannot make " + out_path + " directory.")

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------
def norm_amplitudes(amplitudes):
    """Make the amplitude values span interval [0,1] better
    """
    norm_amps = 1 + 0.1 * (np.log10(np.array(amplitudes)))
    return [max([x, 0]) for x in norm_amps]
    #amplitudes = np.array(amplitudes)
    #return amplitudes / max(amplitudes)

#=============================================================================
# Loop through tests
#=============================================================================
for itest, test in enumerate(tests):
    ntest = itest + 1
    print('Test ' + str(ntest))

    #=========================================================================
    # Models for analysis
    #=========================================================================
    if ntest == 1:
        #---------------------------------------------------------------------
        # Test 1. effect of odor vs. no odor (asleep, maximum concentration)
        #---------------------------------------------------------------------
        rows = [8, 9]
        conditions = [0, 0]
        amplitudes = [1, 1]
        onsets = [onset_list[0], onset_list[1]]
        durations = duration_list
    elif ntest == 2:
        #---------------------------------------------------------------------
        # Test 2. effect of odor vs. no odor (awake, maximum concentration)
        #---------------------------------------------------------------------
        rows = [18, 19]
        conditions = [0, 0]
        amplitudes = [1, 1]
        onsets = [onset_list[0], onset_list[1]]
        durations = duration_list
    elif ntest == 3:
        #---------------------------------------------------------------------
        # Test 3. effect of concentration (asleep)
        #---------------------------------------------------------------------
        rows = range(2, 10)
        n_runs = len(rows)
        conditions = np.zeros(2 * len(amplitude_list), dtype=int).tolist()
        conditions.extend([x + 1 for x in range(n_runs)])
        oneruns = [1 for x in range(n_runs)]
        amplitudes = [[x, x] for x in amplitude_list]
        amplitudes.append(oneruns)
        amplitudes = [x for lst in amplitudes for x in lst]
        # Normalize amplitudes
        amplitudes = norm_amplitudes(amplitudes)
        durations = []
        [durations.extend(duration_list) for x in range(n_runs)]
        onsets = []
        for irun in range(n_runs):
            offset = irun * images_per_run
            onsets.append(offset + onset_list[0])
            onsets.append(offset + onset_list[1])
    elif ntest == 4:
        #---------------------------------------------------------------------
        # Test 4. effect of concentration (awake)
        #---------------------------------------------------------------------
        rows = range(12, 20)
        n_runs = len(rows)
        conditions = np.zeros(2 * len(amplitude_list), dtype=int).tolist()
        conditions.extend([x + 1 for x in range(n_runs)])
        oneruns = [1 for x in range(n_runs)]
        amplitudes = [[x, x] for x in amplitude_list]
        amplitudes.append(oneruns)
        amplitudes = [x for lst in amplitudes for x in lst]
        # Normalize amplitudes
        amplitudes = norm_amplitudes(amplitudes)
        durations = []
        [durations.extend(duration_list) for x in range(n_runs)]
        onsets = []
        for irun in range(n_runs):
            offset = irun * images_per_run
            onsets.append(offset + onset_list[0])
            onsets.append(offset + onset_list[1])
    elif ntest == 5:
        #---------------------------------------------------------------------
        # Test 5. effect of asleep vs. awake (all concentrations)
        #---------------------------------------------------------------------
        rows_asleep = range(2, 10)
        rows_awake = range(12, 20)
        rows = rows_asleep
        rows.extend(rows_awake)
        conditions = [0, 1]
        amplitudes = [1, 1]
        onsets = [0, len(rows_asleep)/2 * images_per_run]
        durations = [len(rows_asleep)/2 * images_per_run,
                     len(rows_awake)/2 * images_per_run]

    #=========================================================================
    # Preprocess (coregister, divide, and smooth) images
    #=========================================================================
    n_images = len(rows)/2 * images_per_run
    converted_file1 = os.path.join(out_path,
        'converted' + wavelengths[0] + '_test' + str(ntest) + ext)
    converted_file2 = os.path.join(out_path,
        'converted' + wavelengths[1] + '_test' + str(ntest) + ext)
    ratio_file = os.path.join(out_path, 'ratio_test' + str(ntest) + ext)
    moco_file =  os.path.join(out_path, 'moco_test' + str(ntest) + ext)
    smooth_file = os.path.join(out_path, 'smooth_test' + str(ntest) + ext)
    #-------------------------------------------------------------------------
    # Convert each .pst image file listed in each row of the table
    # to nifti (neuroimaging file) format and create a slice stack
    # corresponding to each of two wavelengths
    #-------------------------------------------------------------------------
    if convert_images:
        print('Convert images...')

        # Load table
        try:
            csv_reader = csv.reader(open(table_file, 'rU'), dialect=csv.excel_tab)
        except IOError:
            print("  Cannot open " + table_file + ".")

        # Loop through rows
        count1 = 0
        count2 = 0
        image_stack1 = np.zeros((xdim, ydim, 1, n_images))
        image_stack2 = np.zeros((xdim, ydim, 1, n_images))
        for irow, row in enumerate(csv_reader):
            if irow in rows:

                # Load .pst file
                file = os.path.join(images_dir, row[image_file_column])
                wavelength = row[wavelength_column]
                print('  Loading ' + file + ' and converting images...')
                raw = np.fromfile(file, dtype='<i2')
                for iframe in range(images_per_run):
                    image_vector = raw[iframe * xdim * ydim : (iframe + 1) * xdim * ydim]
                    image_matrix = np.reshape(image_vector, (xdim, ydim))

                    # Stack images
                    if wavelength == wavelengths[0]:
                        image_stack1[:, :, 0, count1] = image_matrix
                        count1 += 1
                    elif wavelength == wavelengths[1]:
                        image_stack2[:, :, 0, count2] = image_matrix
                        count2 += 1

        nb.save(nb.Nifti1Image(image_stack1, np.eye(4)), converted_file1)
        nb.save(nb.Nifti1Image(image_stack2, np.eye(4)), converted_file2)

    #-------------------------------------------------------------------------
    # Divide the image volume for one wavelength by the motion-corrected
    # image for the second wavelength (assumded to be aligned slice-wise)
    #-------------------------------------------------------------------------
    if divide_images:
        print('Dividing image volume of wavelength 1 by image volume of wavelength 2...')
        cmd = ['  fslmaths', converted_file1, '-div', converted_file2, ratio_file]
        print(' '.join(cmd)); os.system(' '.join(cmd))

    #-------------------------------------------------------------------------
    # Apply FSL's motion correction to each slice stack
    #-------------------------------------------------------------------------
    if correct_motion:
        print('Correcting motion...')
        cmd = ['  mcflirt -in', ratio_file, '-out', moco_file]
        print(' '.join(cmd)); os.system(' '.join(cmd))
        # cmd = ['  mcflirt -in', ratio_file, '-out moco -edge']
        # print(' '.join(cmd)); os.system(' '.join(cmd))
        # cmd = ['  mv crefvol_moco.nii.gz', moco_file]
        # print(' '.join(cmd)); os.system(' '.join(cmd))

    #-------------------------------------------------------------------------
    # Smooth each ratio image with a Gaussian kernel.
    #-------------------------------------------------------------------------
    if smooth_images:
        print('Smoothing image')
        cmd = ['  fslmaths', moco_file, '-s', str(smooth_sigma), smooth_file]
        print(' '.join(cmd)); os.system(' '.join(cmd))

    #=========================================================================
    # Conduct a general linear model analysis on the preprocessed images per test
    # (Requires image_stack and the following paradigm lists from above:
    #  conditions, onsets, durations, amplitudes)
    #=========================================================================
    if run_analysis:
        ('Run general linear model analysis for each test...')
        if not smooth_images:
            img = nb.load(smooth_file)
            image_stack = img.get_data()

        #-----------------------------------------------------------------
        # Construct a design matrix for each test
        #-----------------------------------------------------------------
        print('  Make design matrix...')
        print('    Conditions:\n      {}'.format(conditions))
        print('    Amplitudes:\n      {}'.format(amplitudes))
        print('    Onsets:\n      {}'.format(onsets))
        print('    Durations:\n      {}'.format(durations))
        paradigm = BlockParadigm(con_id=conditions, onset=onsets,
                                 duration=durations, amplitude=amplitudes)
        frametimes = np.linspace(0, n_images-1, n_images)

        dmtx = make_dmtx(frametimes, paradigm, hrf_model='FIR',
                         drift_model='Blank', hfcut=np.inf)
        design_matrix = dmtx.matrix

        # Plot the design matrix
        if plot_design_matrix:
            fig1 = mp.figure(figsize=(10, 6))
            dmtx.show()
            mp.title('Block design matrix')
            fig1_file = os.path.join(out_path, 'design_matrix_test' + str(ntest) + '.png')
            mp.savefig(fig1_file)
        """
        #-----------------------------------------------------------------
        # Apply a general linear model to all pixels
        #-----------------------------------------------------------------
        print('   Apply general linear model...')
        method = "kalman"
        model = "ar1"
        glm = GLM.glm()
        glm.fit(image_stack.T, design_matrix, method=method, model=model)

        #-----------------------------------------------------------------
        # Create a contrast image
        #
        # Contrast condition 1 vs. condition 2, holding condition 3 constant
        # (sleep vs. awake holding concentration of odorant constant)
        #-----------------------------------------------------------------
        print('  Make contrast image...')

        # Specify the contrast [1 -1 0 ..]
        contrast = np.zeros(design_matrix.shape[1])
        contrast[0] = 1
        contrast[1] = -1
        glm_contrast = glm.contrast(contrast)

        # Compute the contrast image
        zvalues = glm_contrast.zscore().squeeze()
        contrast_image = nb.Nifti1Image(zvalues.T, np.eye(4))

        # Save the contrast as an image
        contrast_file = os.path.join(out_path, 'zmap_test' + str(ntest) + ext)
        nb.save(contrast_image, contrast_file)

        # Plot histogram
        if plot_histogram:
            h, c = np.histogram(zvalues, 100)
            fig2 = mp.figure()
            mp.plot(c[: - 1], h)
            mp.title('Test' + str(ntest) + ': Histogram of the z-values')
            fig2_file = os.path.join(out_path, 'histogram_test' + str(ntest) + '.png')
            mp.savefig(fig2_file)

        # Plot contrast image
        if plot_contrast:
            fig3 = mp.figure()
            mp.matshow(zvalues)
            mp.title('Test' + str(ntest) + ': Contrast image')
            fig3_file = os.path.join(out_path, 'contrast_test' + str(ntest) + '.png')
            mp.savefig(fig3_file)
        """