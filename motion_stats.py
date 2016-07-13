#!/usr/bin/env python3
"""
Extract motion stats from a FSL FEAT or ICA directory

Usage
----
motion_stats.py -i <FSL motion parameter file>
motion_stats.py -h

Authors
----
Mike Tyszka, Caltech, Division of Humaninities and Social Sciences

Dates
----
2016-05-22 JMT Adapt from old Matlab code

License
----
This file is part of MRIutils.

    MRIutils is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRIutils is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRIutils.  If not, see <http://www.gnu.org/licenses/>.

Copyright
----
2016 California Institute of Technology.
"""

__version__ = '0.1.0'

import os
import sys
import argparse
import numpy as np


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Summarize motion stats from preprocessed FEAT directory')
    parser.add_argument('-i', '--mcf_file', required=True, help="FSL MCFLIRT motion parameter file (*_mcf.par)")

    # Parse command line arguments
    args = parser.parse_args()
    mcf_file = args.mcf_file

    # Check for existance
    if not os.path.isfile(mcf_file):
        print('*** FSL MCFLIRT motion parameter file not found (%s)' % mcf_file)
        sys.exit(1)

    # Read file into 2D array
    try:
        M = np.loadtxt(mcf_file)
    except:
        print('*** Could not load motion data from file')
        sys.exit(1)

    # Extract cumulative rotations and displacements
    rx, ry, rz = M[:,0], M[:,1], M[:,2]
    x, y, z = M[:,3], M[:,4], M[:,5]

    # Frame-to-frame displacements and rotations
    drx, dry, drz = np.diff(rx), np.diff(ry), np.diff(rz)
    dx, dy, dz = np.diff(x), np.diff(y), np.diff(z)

    # Total F-F rotations and displacements
    dr, _ = total_rotation(drx, dry, drz)
    dd = np.sqrt(dx**2 + dy**2 + dz**2)

    # Mean F-F rotations (mdeg) and displacements (microns)
    mean_dr = np.mean(dr) * 1000.0 * 180.0 / np.pi
    mean_dd = np.mean(dd) * 1000.0

    # Output mean F-F rotation (mdeg) and displacement (microns)
    print(mean_dr, mean_dd)


def total_rotation(rx, ry, rz):
    """
    Total rotation and axes for a three-axis gimble rotation
    Adapted from the discussion thread at http://www.mathworks.com/matlabcentral/newsreader/view_thread/160945

    :param rx: numpy array
    :param ry: numpy array
    :param rz: numpy array
    :return: rtot, raxes

    """
    # Find length of angle vectors
    n = rx.shape[0]

    # Make space for results
    rtot = np.zeros([n, 1])
    raxes = np.zeros([n,3])

    for i, tx in enumerate(rx):

        # Construct total rotation matrix
        ctx = np.cos(tx)
        stx = np.sin(tx)
        Rx = np.array([[1,0,0], [0,ctx,-stx], [0,stx,ctx]])

        cty = np.cos(ry[i])
        sty = np.sin(ry[i])
        Ry = np.array([[cty,0,sty], [0,1,0] ,[-sty,0,cty]])

        ctz = np.cos(rz[i])
        stz = np.sin(rz[i])
        Rz = np.array([[ctz,-stz,0], [stz,ctz,0], [0,0,1]])

        # Total rotation matrix is product of axis rotations
        Rtot = np.dot(Rz, np.dot(Ry,Rx))

        # Direct calculation of angle and axis from A
        # Code adapted from thread response by Bruno Luong

        # Rotation axis u = [x, y, z]
        u = np.array([Rtot[2,1]-Rtot[1,2], Rtot[0,2]-Rtot[2,0], Rtot[1,0]-Rtot[0,1]])

        # Rotation sine and cosine
        c = np.trace(Rtot) - 1
        s = np.linalg.norm(u)

        # Total rotation in radians
        rtot[i] = np.arctan2(s, c)

        # Adjust rotation to be positive, flipping axis if necessary
        if s > 0:
            u /= s
        else:
            # warning('A close to identity, arbitrary result');
            u = [1, 0, 0]

        # Save axis result
        raxes[i, :] = u

    return rtot, raxes


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()