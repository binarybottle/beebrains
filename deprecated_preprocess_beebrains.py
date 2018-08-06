#!/usr/bin/python

"""
(1) Open each .lst (csv) table in a given directory.
(2) Open .pst image files listed in each row of the table.
(3) Divide each .pst image corresponding to one 
    wavelength by the image of a second wavelength.
(4) To motion-correct the images from (3), register each 
    image to the middle image by computing an affine then
    a nonlinear transform (ANTS software from UPenn).
(5) Smooth each motion-corrected image with a sigma=3 Gaussian kernel.
(6) Save individual motion-corrected and smoothed+motion-corrected images,
    as well as concatenated images from (4) and from (5),
    as nifti files (.nii.gz) for use with fMRI brain imaging software.

Outputs: Nifti files for each table (for each bee).

Command: python <this file name> <data directory> <output directory>

Example: python preprocess_beebrains.py data output

Requirements:  
* Python libraries:  nibabel, numpy, scipy
* ANTS registration software for motion correction
* ImageMagick for converting images and creating montage movies

(c) 2012  Arno Klein, under Apache License Version 2.0
          arno@binarybottle.com  .  www.binarybottle.com
"""

import os, sys
from glob import glob
import csv
import nibabel as nib
import numpy as np
from scipy.ndimage.filters import gaussian_filter

# Settings
xdim = 130
ydim = 172
frames_per_run = 232
#frame_stim1_on = '72'
#frame_stim1_off = '84'
#frame_stim2_on = '92'
#frame_stim2_off = '104'

# Steps to run:
convert_images = 1
correct_motion = 1
smooth_images = 1

# Save image types:
concatenate_images = 1
save_each_image = 1
save_affine = 1
save_nonlinear = 1
# Save for visualizing results:
save_montages = 0
save_movie = 0

# Command-line arguments
args = sys.argv[1:]
if len(args)<2:
    print "\n\t Please provide the names of two directories: one containing .lst table files, another to save output."
    print "\t Example: python " + sys.argv[0] + " data output"
    sys.exit()

# Loop through table files
in_path = str(args[0])
out_path = str(args[1])
try:
    os.system('mkdir ' + out_path)
except IOError:
    print "Cannot make " + out_path + " directory."
try:
    table_files = glob(os.path.join(in_path,"*.lst"))
except IOError:
    print "Cannot make sense of provided path."
for itable_file in table_files:
    table_file = str(itable_file)
    if os.path.exists(table_file):
        print('Table file:  ' + table_file)
        output_stem = os.path.join(out_path, table_file.split('/')[-1])
        output_path_images = output_stem + '_images'

        # Load table and count how many rows are to be included (analysis = 2)
        spamReader = csv.reader(open(table_file, 'r'), delimiter="\t")
        count_rows = 0
        for row in spamReader:
            if row[-1] == '2':
                count_rows += 1

        frames_total = count_rows*frames_per_run

        # Convert images to nifti files
        if convert_images:
            if concatenate_images:
                print('Initializing volume for concatenating image files...')
                B = np.zeros((xdim, ydim, frames_total))            
            if save_each_image:
                os.system('mkdir ' + output_path_images)

            # Load table
            spamReader = csv.reader(open(table_file, 'r'), delimiter="\t")

            # Loop through rows
            count = 0
            for row in spamReader:
                if row[-1] == '2':
                
                    # Extract elements from row
                    file_lambda1 = os.path.join(in_path, row[3].replace('\\','/')) + '.pst'
                    file_lambda2 = os.path.join(in_path, row[18].replace('\\','/')) + '.pst'
                    odor = row[2]
                    state = row[17]

                    # Print elements from row
                    print('File wavelength 1:  ' + file_lambda1)
                    print('File wavelength 2:  ' + file_lambda2)
                    #print('Stimulus frames:  ' + \
                    #       frame_stim1_on + '-' + frame_stim1_off + \
                    #       ' and ' + \
                    #       frame_stim2_on + '-' + frame_stim2_off)
                    print('State:  ' + state)
                    print('Odor:  ' + odor)
            
                    # Load .pst files for both wavelengths
                    print('Loading .pst files...')
                    raw1 = np.fromfile(file_lambda1, dtype='<i2')
                    raw2 = np.fromfile(file_lambda2, dtype='<i2')

                    # Divide the first image by the second
                    print('Dividing images for each of two wavelengths...')
                    raw = (1.0 * raw1) / raw2

                    # Loop through images
                    print('Concatenating resulting images...')
                    for iframe in range(frames_per_run):
                        image_vector = raw[iframe*xdim*ydim : (iframe+1)*xdim*ydim]
                        image_matrix = np.reshape(image_vector, (xdim,ydim))
    
                        # Concatenate images to create one matrix per table
                        if concatenate_images:
                            B[:, :, count] = image_matrix
 
                        count += 1

                        # Save each ratio image as a nifti image
                        if save_each_image:
                            img_ratio_nib = nib.Nifti1Image(image_matrix, np.eye(4))
                            img_ratio_nib.to_filename(os.path.join(output_path_images,'image'+str(count)+'.nii.gz'))
  
            # Save the concatenated nifti file
            print('Saving concatenated image volume')
            if concatenate_images:
                B_nib = nib.Nifti1Image(B, np.eye(4))
                B_nib.to_filename(output_stem + '.nii.gz')

        # Apply motion correction to first image
        # Warp: ANTS 3 -m CC[target.nii.gz, source.nii.gz, 1, 2] -o output_transform.nii.gz 
        #              -r Gauss[2,0] -t SyN[0.5] -i 30x99x11 --use-Histogram-Matching 
        #              --number-of-affine-iterations 10000x10000x10000x10000x10000 
        # Reslice: WarpImageMultiTransform 3 labeled_source.nii.gz output_labels.nii.gz 
        #              -R target.nii transformWarp.nii.gz transformAffine.txt -use-NN 
        if correct_motion:
            print('Correcting motion...')
            iframe_middle = np.int(frames_total/2)
            image_ref = os.path.join(output_path_images,'image'+str(iframe_middle)+'.nii.gz')
            for iframe in range(frames_total):
                image = os.path.join(output_path_images,'image'+str(iframe+1)+'.nii.gz')
                image_mc_stem = os.path.join(output_path_images,'image'+str(iframe+1)+'_motioncorrected')
                image_mc = image_mc_stem + '.nii.gz'
                cmd1 = ['ANTS 2 -m CC[',image_ref+',',image+',1,2] -o',image_mc,
                        '-r Gauss[2,0] -t SyN[0.5] -i 30x99x11 --use-Histogram-Matching',
                        '--number-of-affine-iterations 10000x10000x10000x10000x10000']
                print(' '.join(cmd1))
                os.system(' '.join(cmd1))
                if save_affine:
                    image_mc_affine = image_mc_stem + 'Affine.nii.gz'
                    cmd2 = ['WarpImageMultiTransform 2',image,image_mc_affine,'-R',image_ref,
                            image_mc_stem+'Affine.txt']
                    print(' '.join(cmd2))
                    os.system(' '.join(cmd2))
                if save_nonlinear:
                    image_mc_nonlinear = image_mc
                    cmd3 = ['WarpImageMultiTransform 2',image,image_mc_nonlinear,'-R',image_ref,
                            image_mc_stem+'Warp.nii.gz',
                            image_mc_stem+'Affine.txt']
                    print(' '.join(cmd3))
                    os.system(' '.join(cmd3))

        # Smooth each motion-corrected image file.
        # Concatenate motion-corrected images.
        # Concatenate smoothed, motion-corrected images.
        if smooth_images:
            if concatenate_images:
                print('Initializing volumes for motion-corrected and smoothed image files...')
                Bc = np.zeros((xdim, ydim, frames_total))                      
                if save_affine:
                    Bc_affine_smooth = Bc.copy()
                if save_nonlinear:
                    Bc_nonlinear_smooth = Bc.copy()
            for iframe in range(frames_total):
                print('Smoothing and concatenating image ' + str(iframe+1))
                image_mc_stem = os.path.join(output_path_images,'image'+str(iframe+1)+'_motioncorrected')
                if save_affine:
                    image_mc_affine = image_mc_stem + 'Affine.nii.gz'
                    image_mc_affine_nib = nib.load(image_mc_affine)
                    image_mc_affine = image_mc_affine_nib.get_data()
                    # Smooth each image
                    image_mc_affine_smooth = gaussian_filter(image_mc_affine, sigma=3, order=0)
                if save_nonlinear:
                    image_mc_nonlinear = image_mc_stem + '.nii.gz'
                    image_mc_nib = nib.load(image_mc_nonlinear)
                    image_mc = image_mc_nib.get_data()
                    # Smooth each image
                    image_mc_smooth = gaussian_filter(image_mc, sigma=3, order=0)

                if concatenate_images:
                    Bc[:, :, iframe] = image_mc
                    if save_affine:
                        Bc_affine_smooth[:, :, iframe] = image_mc_affine_smooth 
                    if save_nonlinear:
                        Bc_nonlinear_smooth[:, :, iframe] = image_mc_smooth 
                if save_each_image:
                    if save_affine:
                        image_mc_affine_smooth_nib = nib.Nifti1Image(image_mc_affine_smooth, np.eye(4))
                        image_mc_affine_smooth_file = 'image'+str(iframe+1)+'_motioncorrectedAffine_smoothed.nii.gz'
                        image_mc_affine_smooth_nib.to_filename(os.path.join(output_path_images,image_mc_affine_smooth_file))
                    if save_nonlinear:
                        image_mc_smooth_nib = nib.Nifti1Image(image_mc_smooth, np.eye(4))
                        image_mc_smooth_file = 'image'+str(iframe+1)+'_motioncorrected_smoothed.nii.gz'
                        image_mc_smooth_nib.to_filename(os.path.join(output_path_images,image_mc_smooth_file))

            # Save the concatenated nifti files
            if concatenate_images:
                Bc_nib = nib.Nifti1Image(Bc, np.eye(4))
                Bc_nib.to_filename(output_stem + '_motioncorrected.nii.gz')
                if save_affine:
                    Bc_affine_smooth_nib = nib.Nifti1Image(Bc_affine_smooth, np.eye(4))
                    Bc_affine_smooth_nib.to_filename(output_stem + '_motioncorrectedAffine_smoothed.nii.gz')
                if save_nonlinear:
                    Bc_nonlinear_smooth_nib = nib.Nifti1Image(Bc_nonlinear_smooth, np.eye(4))
                    Bc_nonlinear_smooth_nib.to_filename(output_stem + '_motioncorrected_smoothed.nii.gz')

        # Convert each nifti image file to jpg and to png and save a montage of preprocessed images
        if save_montages:
            for iframe in range(frames_total):
                print('Converting nifti to jpg and png image and creating a montage for' + str(iframe+1))
                image_stem = os.path.join(output_path_images,'image'+str(iframe+1))
                output_montage = os.path.join(output_path_images,'montage'+str(iframe+1))+'.png'
                cmd = ['ConvertToJpg',image_stem+'.nii.gz',image_stem+'.jpg']
                print(' '.join(cmd)); os.system(' '.join(cmd))                
                cmd = ['convert',image_stem+'.jpg',image_stem+'.png']
                print(' '.join(cmd)); os.system(' '.join(cmd))                
                image_mc_stem = os.path.join(output_path_images,'image'+str(iframe+1)+'_motioncorrected')
                if save_affine:
                    image_mc_affine_stem = image_mc_stem + 'Affine'
                    cmd = ['ConvertToJpg',image_mc_affine_stem+'.nii.gz',image_mc_affine_stem+'.jpg']
                    print(' '.join(cmd)); os.system(' '.join(cmd))                
                    cmd = ['convert',image_mc_affine_stem+'.jpg',image_mc_affine_stem+'.png']
                    print(' '.join(cmd)); os.system(' '.join(cmd))                
                    image_mc_affine_sm_stem = image_mc_affine_stem + '_smoothed'
                    cmd = ['ConvertToJpg',image_mc_affine_sm_stem+'.nii.gz',image_mc_affine_sm_stem+'.jpg']
                    print(' '.join(cmd)); os.system(' '.join(cmd))                
                    cmd = ['convert',image_mc_affine_sm_stem+'.jpg',image_mc_affine_sm_stem+'.png']
                    print(' '.join(cmd)); os.system(' '.join(cmd))                
                if save_nonlinear:
                    cmd = ['ConvertToJpg',image_mc_stem+'.nii.gz',image_mc_stem+'.jpg']
                    print(' '.join(cmd)); os.system(' '.join(cmd))                
                    cmd = ['convert',image_mc_stem+'.jpg',image_mc_stem+'.png']
                    print(' '.join(cmd)); os.system(' '.join(cmd))                
                    image_mc_sm_stem = image_mc_stem + '_smoothed'
                    cmd = ['ConvertToJpg',image_mc_sm_stem+'.nii.gz',image_mc_sm_stem+'.jpg']
                    print(' '.join(cmd)); os.system(' '.join(cmd))                
                    cmd = ['convert',image_mc_sm_stem+'.jpg',image_mc_sm_stem+'.png']
                    print(' '.join(cmd)); os.system(' '.join(cmd))                
                cmd = ['montage',
                       image_stem+'.png',
                       image_mc_affine_stem+'.png',
                       image_mc_stem+'.png',
                       image_mc_sm_stem+'.png',
                       '-label ' + str(iframe+1),
                       '-geometry +1+0 -tile 4x -background black',
                       output_montage]
                print(' '.join(cmd)); os.system(' '.join(cmd))                
        if save_movie:
            input_montages = os.path.join(output_path_images,'montage*.png')
            output_montage_movie = output_stem + '_montages.gif'
            cmd = ['convert',input_montages,'-adjoin -compress none',output_montage_movie]
            print(' '.join(cmd)); os.system(' '.join(cmd))                
