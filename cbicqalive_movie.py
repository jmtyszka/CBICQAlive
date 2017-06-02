#!/usr/bin/env python3
#
# Create a movie of orthoslices through a 4D fMRI dataset
#
# USAGE : cbicqalive_movie.py -i <4D fMRI> -o <MP4 movie> -c <slice intersection in form x,y,z>
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2017-05-19 JMT From scratch
#
# This file is part of CBICQAlive.
#
#    CBICQAlive is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    CBICQAlive is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#   along with CBICQAlive.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2017 California Institute of Technology.

import os
import sys
import argparse
import nibabel as nib
import numpy as np
import cv2
import datetime
from skimage.exposure import rescale_intensity
from skimage.transform import warp, SimilarityTransform


# Main function
def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='CBICQAlive statistics')
    parser.add_argument('-i', '--fmri', required=True, help="4D compressed Nifti fMRI dataset")
    parser.add_argument('-c', '--center', nargs="+", type=int, help="Intersection coordinates eg 57 48 19 [3D image center]")

    # Parse command line arguments
    args = parser.parse_args()

    infile = args.fmri

    # Load the fMRI dataset
    try:
        print('  Loading image from %s' % infile)
        nii = nib.load(infile)
        img = nii.get_data()
    except:
        print('* Problem loading data from %s' % infile)
        print('* Exiting')
        sys.exit(1)

    # Init the output directory
    out_dir = infile.replace('.nii.gz', '.tryptic')
    print('  All tryptics will be saved in %s' % out_dir)
    os.makedirs(out_dir, exist_ok=True)

    # fMRI data dimensions
    nx, ny, nz, nt = img.shape

    # TR in seconds from header
    TR_s = nii.header.get_zooms()[3]

    if args.center:
        x0, y0, z0 = args.center
    else:
        x0, y0, z0 = int(nx/2), int(ny/2), int(nz/2)

    print('  Intersection set at (%d, %d, %d)' % (x0, y0, z0))

    # Set timestamp font and color
    tstamp_font = cv2.FONT_HERSHEY_SIMPLEX
    tstamp_color = (255, 255, 255)

    # Use global 99th percentile as upper limit for robust rescaling
    robust_range = (0, np.int(np.percentile(img, 99.0)))
    print('  Robust rescaling to range [%d, %d]' % robust_range)

    # Tryptic height in pixels (k)
    # Tryptic width = 3k
    k = 400

    # Voxel scale factor based on maximum spatial dimension
    sf = float(k) / float(np.max([nx, ny, nz]))
    print('  Voxels scaled by %0.3f' % sf)

    # Loop over all time points
    for tc in range(0,nt):

        print('  Processing volume %04d' % tc)

        # Extract xy, xz and yz slices through intersection
        img_xy = np.rot90(img[:,:,z0,tc])
        img_xz = np.rot90(img[:,y0,:,tc])
        img_yz = np.rot90(img[x0,:,:,tc])

        # Resize slices to constant height, holding aspect ratio constant
        img_xy = myresize(img_xy, sf=sf, k=k)
        img_xz = myresize(img_xz, sf=sf, k=k)
        img_yz = myresize(img_yz, sf=sf, k=k)

        # Create tryptic and robust rescale to unit8 full range [0,255]
        tryptic = np.hstack([img_xy, img_xz, img_yz])
        tryptic = np.uint8(rescale_intensity(tryptic, in_range=robust_range, out_range=np.uint8))

        # Add colored timestamp text
        tryptic_rgb = cv2.cvtColor(tryptic, cv2.COLOR_GRAY2RGB)
        tstamp_str = str(datetime.timedelta(seconds=tc * TR_s))
        cv2.putText(tryptic_rgb, tstamp_str, (10, 40), tstamp_font, 1, tstamp_color, 2)

        cv2.imshow("Tryptic", tryptic_rgb)

        # Save image
        out_fname = os.path.join(out_dir, 'tryptic_%04d.png' % tc)
        cv2.imwrite(out_fname, tryptic_rgb)

        cv2.waitKey(1)

    sys.exit(0)


def myresize(img, sf=1.0, k=640):

    # Center scaled image in output frame
    hk, hx, hy = k/2.0, img.shape[1]/2.0, img.shape[0]/2.0
    dx = hx - hk/sf
    dy = hy - hk/sf

    tx = SimilarityTransform(scale=1/sf, translation=(dx,dy))
    img_resize = warp(img, tx, output_shape=(k, k), order=0, preserve_range=True)

    return img_resize


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
