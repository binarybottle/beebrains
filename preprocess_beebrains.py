#!/usr/bin/python

"""
Preprocessing steps:

(1) Open each .lst (csv) table in a given directory.
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
(6) Smooth each ratio image (from 5) with a sigma=3 Gaussian kernel.

Outputs: Nifti files for each table (for each bee).

Command: python <this file name> <data directory> <output directory>

Example: python preprocess_beebrains.py data output

Requirements:  
* Python libraries:  nibabel, numpy, scipy
* ANTS registration software for motion correction
* ImageMagick -- only if converting images and creating montage movies

(c) 2012  Arno Klein, under Apache License Version 2.0
          arno@binarybottle.com  .  www.binarybottle.com
"""

import os, sys
from glob import glob
import csv
import nibabel as nib
import numpy as np
from scipy.ndimage.filters import gaussian_filter

#ANTS = '/software/ANTS_1.9/bin/'  # for registration (motion correction)
#IMAGEMAGICK = '/usr/local/bin/'  # only if converting images and creating montage movies
#ANTS = '/software/ANTS/bin/'  # for registration (motion correction)
ANTS = '/usr/local/bin/'  # for registration (motion correction)
IMAGEMAGICK = '/usr/local/bin/'  # only if converting images and creating montage movies

# Settings
xdim = 130  # x dimension for each image
ydim = 172  # y dimension for each image
frames_per_run = 232  # number of images captured for each run

# Steps to run:
convert_images = 0  # convert .pst images to nifti file format
<<<<<<< HEAD
correct_motion = 0  # apply registration to correct for motion
smooth_images = 0  # smooth the resulting motion-corrected images

if smooth_images:
    smooth_sigma = 1

# Table parameters:
include_string = '2'  # indicator in table which rows are to be processed
table_file_string = '*.lst'  # string identifying table files
table_delimiter = "\t"  # the delimiter between elements of the table
#odor_column = 2
#state_column = 17
lambda_columns = [3,18]  # table columns that give the file names for the lambda 1 and 2 images
wavelengths = ['340nm','380nm']  # lambda 1 and 2 strings

#------------------------------------------------------------------------------------------------

# Save image types:
save_affine = 0  #  Save affine transform -- NOTE: Using average of affine transforms instead
save_nonlinear = 0  #  Save nonlinear transform -- NOTE: Using average of nonlinear transforms instead
save_affine_avg = 1  # Use average of affine transforms instead
save_nonlinear_avg = 1  # Use average of nonlinear transforms instead
save_montages = 0  # Save for visualizing results:
save_movie = 0  # Very memory intensive -- better to use another program such as GraphicConverter

def convert2png(stem):
    '''Convert nifti file to jpeg and png'''
    cmd = [ANTS+'ConvertToJpg', stem + '.nii.gz', stem + '.jpg']
    print(' '.join(cmd)); os.system(' '.join(cmd))                
    cmd = [IMAGEMAGICK+'convert', stem + '.jpg', stem + '.png']
    print(' '.join(cmd)); os.system(' '.join(cmd))                

def smooth2d(image_in, image_out, smooth_sigma=smooth_sigma):
    '''Smooth each image with a Gaussian filter, sigma=3'''
    image_nib = nib.load(image_in)
    image = image_nib.get_data()
    image_smooth = gaussian_filter(image, sigma=smooth_sigma, order=0)
    image_smooth_nib = nib.Nifti1Image(image_smooth, np.eye(4))
    image_smooth_nib.to_filename(image_out)

w1 = wavelengths[0]
w2 = wavelengths[1]

# Command-line arguments
args = sys.argv[1:]
if len(args)<2:
    print("\n\t Please provide the names of two directories: one containing .lst table files, another to save output.")
    print("\t Example: python " + sys.argv[0] + " data output")
    sys.exit()

# Loop through table files
in_path = str(args[0])
out_path = str(args[1])
try:
    if not os.path.exists(out_path):
        os.mkdir(out_path)
except IOError:
    print("Cannot make " + out_path + " directory.")
try:
    table_files = glob(os.path.join(in_path, table_file_string))
except IOError:
    print("Cannot make sense of provided path.")
for itable_file in table_files:
    table_file = str(itable_file)
    if os.path.exists(table_file):
        print('Table file:  ' + table_file)
        output_stem = os.path.join(out_path, table_file.split('/')[-1])
        output_path_images = output_stem + '_images'
        output_path_transforms = output_stem + '_transforms'
        output_path_transformed = output_stem + '_transformed'
        output_path_smoothed = output_stem + '_smoothed' + str(smooth_sigma)
        output_path_montages = output_stem + '_montages'

        # Load table and  how many rows are to be included (analysis = 2)
        try:
            spamReader = csv.reader(open(table_file, 'r'), delimiter=table_delimiter)
        except IOError:
            print("Cannot open " + table_file + ".")
        count_rows = 0
        for row in spamReader:
            if row[-1] == include_string:
                count_rows += 1

        frames_total = count_rows*frames_per_run

        # Convert images to nifti files
        if convert_images:
            if not os.path.exists(output_path_images):
                os.mkdir(output_path_images)

            # Loop through lambdas
            for ilambda in range(len(lambda_columns)):
                count = 0
                # Load table
                spamReader = csv.reader(open(table_file, 'r'), delimiter=table_delimiter)
                # Loop through rows
                for row in spamReader:
                    if row[-1] == include_string:
                        #odor = row[odor_column]
                        #state = row[state_column]

                        # Load .pst file
                        file_lambda = os.path.join(in_path, row[lambda_columns[ilambda]].replace('\\','/')) + '.pst'
                        print('Loading ' + file_lambda + ' and converting images...')
                        raw = np.fromfile(file_lambda, dtype='<i2')

                        # Loop through images
                        for iframe in range(frames_per_run):
                            image_vector = raw[iframe*xdim*ydim : (iframe+1)*xdim*ydim]
                            image_matrix = np.reshape(image_vector, (xdim,ydim))

                            count += 1

                            # Save each image as a nifti file
                            converted_file = 'image' + str(count) + '_' + wavelengths[ilambda] + '.nii.gz'
                            img_ratio_nib = nib.Nifti1Image(image_matrix, np.eye(4))
                            img_ratio_nib.to_filename(os.path.join(output_path_images, converted_file))

        # Apply motion correction (ANTS) to middle image
        # Warp: ANTS 2 -m CC[target.nii.gz, source.nii.gz, 1, 2] -o transform.nii.gz 
        #              -r Gauss[2,0] -t SyN[0.5] -i 30x99x11 --use-Histogram-Matching 
        #              --number-of-affine-iterations 10000x10000x10000x10000x10000 
        # Reslice: WarpImageMultiTransform 2 source.nii.gz source2target.nii.gz 
        #              -R target.nii transformWarp.nii.gz transformAffine.txt
        if correct_motion:
            print('Correcting motion...')
            if not os.path.exists(output_path_transforms):
                os.mkdir(output_path_transforms)
            if not os.path.exists(output_path_transformed):
                os.mkdir(output_path_transformed)
            iframe_middle = np.int(frames_total/2)
            for iframe in range(frames_total):
                image_stem = os.path.join(output_path_images, 'image' + str(iframe+1))
                image_ref_stem = os.path.join(output_path_images, 'image' + str(iframe_middle+1))
                transform_stem = os.path.join(output_path_transforms, 'image' + str(iframe+1))
                transformed_stem = os.path.join(output_path_transformed, 'image' + str(iframe+1))
                compute_transform = 1
                if compute_transform:
                    for ilambda in range(len(lambda_columns)):
                        w = '_' + wavelengths[ilambda]
                        image_ref = image_ref_stem + w + '.nii.gz'
                        image_lambda = image_stem + w + '.nii.gz'
                        xfm = transform_stem + w + '_.nii.gz'
                        cmd = [ANTS+'ANTS 2 -m CC[',image_ref+',',image_lambda+',1,2] -o',xfm,
                               '-r Gauss[2,0] -t SyN[0.5] -i 30x99x11 --use-Histogram-Matching',
                               '--number-of-affine-iterations 10000x10000x10000x10000x10000']
                        print(' '.join(cmd)); os.system(' '.join(cmd))

                        if save_affine:
                            cmd = [ANTS+'WarpImageMultiTransform 2', image_lambda, transformed_stem + w + '_Affined.nii.gz',
                                   '-R', image_ref, transform_stem + w + '_Affine.txt']
                            print(' '.join(cmd)); os.system(' '.join(cmd))

                        if save_nonlinear:
                            cmd = [ANTS+'WarpImageMultiTransform 2', image_lambda, transformed_stem + w + '_Warped.nii.gz',
                                   '-R', image_ref, transform_stem + w + '_Warp.nii.gz', transform_stem + w + '_Affine.txt']
                            print(' '.join(cmd)); os.system(' '.join(cmd))

                if save_affine:
                    # Divide the first affine motion-corrected image by the second
                    print('Dividing affine motion-corrected images for each of two wavelengths...')
                    cmd = [ANTS+'ImageMath 2',transformed_stem + '_Affined_ratio.nii.gz','/',
                           transformed_stem + '_' + w1 + '_Affined.nii.gz',
                           transformed_stem + '_' + w2 + '_Affined.nii.gz']
                    print(' '.join(cmd)); os.system(' '.join(cmd))

                if save_nonlinear:
                    # Divide the first affine motion-corrected image by the second
                    print('Dividing nonlinear motion-corrected images for each of two wavelengths...')
                    cmd = [ANTS+'ImageMath 2',transformed_stem + '_Warped_ratio.nii.gz','/',
                           transformed_stem + '_' + w1 + '_Warped.nii.gz',
                           transformed_stem + '_' + w2 + '_Warped.nii.gz']
                    print(' '.join(cmd)); os.system(' '.join(cmd))

                # Compute the average of the lambda files' two affine transforms,
                # then apply it to the ratio of the two lambda files
                if save_affine_avg or save_nonlinear_avg:
                    cmd = [ANTS+'AverageAffineTransform 2',
                           transform_stem + '_AvgAffine.txt',
                           transform_stem + '_' + w1 + '_Affine.txt',
                           transform_stem + '_' + w2 + '_Affine.txt']
                    print(' '.join(cmd)); os.system(' '.join(cmd))
                if save_affine_avg:
                    for ilambda in range(len(lambda_columns)):
                        w = '_' + wavelengths[ilambda]
                        cmd = [ANTS+'WarpImageMultiTransform 2',
                               image_stem + w + '.nii.gz',
                               transformed_stem + w + '_AvgAffined.nii.gz', '-R',
                               image_ref_stem + w + '.nii.gz',
                               transform_stem + '_AvgAffine.txt']
                        print(' '.join(cmd)); os.system(' '.join(cmd))

                    # Divide the first average-affine motion-corrected image by the second
                    print('Dividing average-affine motion-corrected images for each of two wavelengths...')
                    cmd = [ANTS+'ImageMath 2',transformed_stem + '_AvgAffined_ratio.nii.gz','/',
                           transformed_stem + '_' + w1 + '_AvgAffined.nii.gz',
                           transformed_stem + '_' + w2 + '_AvgAffined.nii.gz']
                    print(' '.join(cmd)); os.system(' '.join(cmd))

                # Compute the average of the lambda files' two nonlinear transforms,
                # then apply them to the ratio of the two lambda files
                if save_nonlinear_avg:
                    cmd = [ANTS+'AverageImages 2', transform_stem + '_AvgWarp.nii.gz', '0',
                                                   transform_stem + '_*_Warp.nii.gz']
                    print(' '.join(cmd)); os.system(' '.join(cmd))

                    for ilambda in range(len(lambda_columns)):
                        w = '_' + wavelengths[ilambda]
                        cmd = [ANTS+'WarpImageMultiTransform 2',
                               image_stem + w + '.nii.gz',
                               transformed_stem + w + '_AvgWarped.nii.gz', '-R',
                               image_ref_stem + w + '.nii.gz',
                               transform_stem + '_AvgWarp.nii.gz',
                               transform_stem + '_AvgAffine.txt']
                        print(' '.join(cmd)); os.system(' '.join(cmd))

                    # Divide the first average-warp motion-corrected image by the second
                    print('Dividing average-warp motion-corrected images for each of two wavelengths...')
                    cmd = [ANTS+'ImageMath 2',transformed_stem + '_AvgWarped_ratio.nii.gz','/',
                           transformed_stem + '_' + w1 + '_AvgWarped.nii.gz',
                           transformed_stem + '_' + w2 + '_AvgWarped.nii.gz']
                    print(' '.join(cmd)); os.system(' '.join(cmd))

        # Smooth each motion-corrected image file.
        if smooth_images:
            if not os.path.exists(output_path_smoothed):
                os.mkdir(output_path_smoothed)
            for iframe in range(frames_total):
                print('Smoothing image ' + str(iframe+1))
                transformed_stem = os.path.join(output_path_transformed, 'image' + str(iframe+1))
                smoothed_stem = os.path.join(output_path_smoothed, 'image' + str(iframe+1))
                #if save_affine:
                #    smooth2d(transformed_stem + '_Affined_ratio.nii.gz',
                #             transformed_stem + '_Affined_ratio_smooth.nii.gz',
                #             smooth_sigma = smooth_sigma)
                #if save_affine_avg:
                #    smooth2d(transformed_stem + '_AvgAffined_ratio.nii.gz',
                #             transformed_stem + '_AvgAffined_ratio_smooth.nii.gz',
                #             smooth_sigma = smooth_sigma)
                #if save_nonlinear:
                #    smooth2d(transformed_stem + '_Warped_ratio.nii.gz',
                #             transformed_stem + '_Warped_ratio_smooth.nii.gz',
                #             smooth_sigma = smooth_sigma)
                if save_nonlinear_avg:
                    smooth2d(transformed_stem + '_AvgWarped_ratio.nii.gz',
                             smoothed_stem + '_AvgWarped_ratio_smooth.nii.gz',
                             smooth_sigma = smooth_sigma)


        # Convert each nifti image file to jpg and to png and save a montage of preprocessed images
        if save_montages:
            if not os.path.exists(output_path_montages):
                os.mkdir(output_path_montages)
            for iframe in range(frames_total):
                print('Converting nifti to jpg and png image and creating a montage for ' + str(iframe+1))
                image_stem = os.path.join(output_path_images, 'image' + str(iframe+1))
                transformed_stem = os.path.join(output_path_transformed, 'image' + str(iframe+1))
                smoothed_stem = os.path.join(output_path_smoothed, 'image' + str(iframe+1))
                output_montage = os.path.join(output_path_montages, 'montage' + str(iframe+1)) + '.png'

                convert2png(image_stem + '_' + w1)
                convert2png(image_stem + '_' + w2)
                if save_affine:
                    convert2png(transformed_stem + '_Affined_ratio')
                if save_nonlinear:
                    convert2png(transformed_stem + '_Warped_ratio')
                    #convert2png(image_stem + '_Warped_ratio_smooth')
                if save_affine_avg:
                    convert2png(transformed_stem + '_AvgAffined_ratio')
                if save_nonlinear_avg:
                    convert2png(transformed_stem + '_AvgWarped_ratio')
                    convert2png(smoothed_stem + '_AvgWarped_ratio_smooth')

                cmd = [IMAGEMAGICK+'montage',
                       image_stem + '_' + w1 + '.png',
                       image_stem + '_' + w2 + '.png',
                       transformed_stem + '_AvgAffined_ratio.png',
                       transformed_stem + '_AvgWarped_ratio.png',
                       smoothed_stem + '_AvgWarped_ratio_smooth.png',
                       '-label ' + str(iframe+1),
                       '-geometry +1+0 -tile 5x -background black',
                       output_montage]
                print(' '.join(cmd)); os.system(' '.join(cmd))

        # Save movie of montages
        # NOTE:  very memory intensive -- better to use another program such as GraphicConverter
        if save_movie:
            input_montages = os.path.join(output_path_montages,'montage*.png')
            output_montage_movie = output_stem + '_montages.gif'
            cmd = [IMAGEMAGICK+'convert', input_montages, '-adjoin -compress none', output_montage_movie]
            print(' '.join(cmd)); os.system(' '.join(cmd))                
