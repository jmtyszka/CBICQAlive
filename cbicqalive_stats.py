#!/usr/bin/env python
#
# Derive various descriptive stats from the preprocessed QA data
# - detrended time series parameters
# - spike counts
# - center of mass
# - apparent motion
#
# USAGE : cbicqa_stats.py <QA Directory>
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 04/01/2013 JMT From scratch
#          10/24/2013 JMT Expand to calculate all stats
#          03/03/2014 JMT Add scanner name argument
#
# This file is part of CBICQA.
#
#    CBICQA is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    CBICQA is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#   along with CBICQA.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2013-2014 California Institute of Technology.

import os
import argparse
from pylab import *


# Main function
def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='CBICQAlive statistics')
    parser.add_argument('-i', '--qa_dir', help="CBICQAlive directory (*.qalive)")

    # Parse command line arguments
    args = parser.parse_args()
    qa_dir = args.qa_dir

    print('  Calculating statistics for ' + qa_dir)
    
    # Load ROI timeseries
    print('    Loading ROI timeseries')
    qa_ts_fname = os.path.join(qa_dir, 'qa_roi_ts.txt')
    if not os.path.isfile(qa_ts_fname):
        print(qa_ts_fname + ' does not exist - exiting')
        sys.exit(0)

    x = np.loadtxt(qa_ts_fname)

    # Parse timeseries into vectors
    brain_t = x[:, 0]
    ghost_t = x[:, 1]
    air_t = x[:, 2]

    # Volume number vector
    nv = len(brain_t)
    v = np.linspace(0, nv-1, nv)

    # Load DVARS timeseries
    print('    Loading DVARS timeseries')
    qa_dvars_fname = os.path.join(qa_dir, 'qa_dvars_ts.txt')
    if not os.path.isfile(qa_dvars_fname):
        print(qa_dvars_fname + ' does not exist - exiting')
        sys.exit(0)

    dvars_t = np.loadtxt(qa_dvars_fname)

    #
    # Analyze motion parameters
    #
    print('    Analyzing motion parameters')
    qa_mcf_parfile = os.path.join(qa_dir, 'qa_mcf.par')

    if not os.path.isfile(qa_mcf_parfile):
        print(qa_mcf_parfile + ' does not exist - exiting')
        sys.exit(0)

    # Read file into 2D array
    try:
        mcf_pars = np.loadtxt(qa_mcf_parfile)
    except:
        print('*** Could not load motion data from file')
        sys.exit(1)

    # Extract cumulative rotations and displacements
    rx, ry, rz = mcf_pars[:, 0], mcf_pars[:, 1], mcf_pars[:, 2]
    x, y, z = mcf_pars[:, 3], mcf_pars[:, 4], mcf_pars[:, 5]

    # Frame-to-frame displacements and rotations
    drx, dry, drz = np.gradient(rx), np.gradient(ry), np.gradient(rz)
    dx, dy, dz = np.gradient(x), np.gradient(y), np.gradient(z)

    # Total F-F rotations and displacements
    dr, _ = total_rotation(drx, dry, drz)
    dd = np.sqrt(dx**2 + dy**2 + dz**2)

    # Convert distances from mm to microns
    dscale = 1000.0
    x_um, y_um, z_um = x * dscale, y * dscale, z * dscale
    dd_um, dx_um, dy_um, dz_um = dd * dscale, dx * dscale, dy * dscale, dz * dscale

    # Convert rotations from radians to millidegrees
    rscale = 180.0 * 1000.0 / np.pi
    rx_mdeg, ry_mdeg, rz_mdeg = rx * rscale, ry * rscale, rz * rscale
    dr_mdeg, drx_mdeg, dry_mdeg, drz_mdeg = dr * rscale, drx * rscale, dry * rscale, drz * rscale

    #
    # Summary statistics for all timeseries
    #
    brain_tmean, brain_thresh, brain_nout, brain_pout, _ = timeseries_stats(brain_t)
    ghost_tmean, ghost_thresh, ghost_nout, ghost_pout, _ = timeseries_stats(ghost_t)
    air_tmean, air_thresh, air_nout, air_pout, _ = timeseries_stats(air_t)
    dvars_tmean, dvars_thresh, dvars_nout, dvars_pout, _ = timeseries_stats(dvars_t)
    dd_um_tmean, dd_um_thresh, dd_um_nout, dd_um_pout, _ = timeseries_stats(dd_um)
    dr_mdeg_tmean, dr_mdeg_thresh, dr_mdeg_nout, dr_mdeg_pout, _ = timeseries_stats(dr_mdeg)

    #
    # Output summary statistics to CSV file
    #
    stats_csv_fname = os.path.join(qa_dir, 'qa_stats.csv')
    fd = open(stats_csv_fname, "w")

    fd.write("Brain, %0.1f, %0.1f, %d, %0.1f\n" % (brain_tmean, brain_thresh, brain_nout, brain_pout))
    fd.write("Ghost, %0.1f, %0.1f, %d, %0.1f\n" % (ghost_tmean, ghost_thresh, ghost_nout, ghost_pout))
    fd.write("Air, %0.1f, %0.1f, %d, %0.1f\n" % (air_tmean, air_thresh, air_nout, air_pout))
    fd.write("DVARS, %0.1f, %0.1f, %d, %0.1f\n" % (dvars_tmean, dvars_thresh, dvars_nout, dvars_pout))
    fd.write("FF_Disp_um, %0.1f, %0.1f, %d, %0.1f\n" % (dd_um_tmean, dd_um_thresh, dd_um_nout, dd_um_pout))
    fd.write("FF_Rot_mdeg, %0.1f, %0.1f, %d, %0.1f\n" % (dr_mdeg_tmean, dr_mdeg_thresh, dr_mdeg_nout, dr_mdeg_pout))

    fd.close()

    # ROI timeseries figure

    fig = figure(figsize=(10, 10))

    subplot(411)
    plot(v, brain_t)
    title("Mean Brain Signal", x=0.5, y=0.8)
    
    subplot(412)
    plot(v, ghost_t)
    title("Mean Ghost Signal", x=0.5, y=0.8)
    
    subplot(413)
    plot(v, air_t)
    title("Mean Air Signal", x=0.5, y=0.8)

    subplot(414)
    plot(v, dvars_t)
    title("DVARS", x=0.5, y=0.8)

    savefig(os.path.join(qa_dir, 'qa_roi_timeseries.png'), dpi=72, bbox_inches='tight')

    # Motion timeseries figure

    fig.clf()

    subplot(411)
    plot(v, x_um, 'r', v, y_um, 'g', v, z_um, 'b')
    title("Displacement (microns)", x=0.5, y=0.8)

    subplot(412)
    plot(v, rx_mdeg, 'r', v, ry_mdeg, 'g', v, rz_mdeg, 'b')
    title("Rotation (mdeg)", x=0.5, y=0.8)

    subplot(413)
    plot(v, dd_um)
    title("Total F-F Displacement (microns)", x=0.5, y=0.8)

    subplot(414)
    plot(v, dr_mdeg)
    title("Total F-F Rotation (mdeg)", x=0.5, y=0.8)

    savefig(os.path.join(qa_dir, 'qa_motion_timeseries.png'), dpi=72, bbox_inches='tight')

    # Done
    print('  Finished python statistical analysis')


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
    raxes = np.zeros([n, 3])

    for i, tx in enumerate(rx):

        # Construct total rotation matrix
        ctx = np.cos(tx)
        stx = np.sin(tx)
        rotx = np.array([[1, 0, 0], [0, ctx, -stx], [0, stx, ctx]])

        cty = np.cos(ry[i])
        sty = np.sin(ry[i])
        roty = np.array([[cty, 0, sty], [0, 1, 0], [-sty, 0, cty]])

        ctz = np.cos(rz[i])
        stz = np.sin(rz[i])
        rotz = np.array([[ctz, -stz, 0], [stz, ctz, 0], [0, 0, 1]])

        # Total rotation matrix is product of axis rotations
        rotall = np.dot(rotz, np.dot(roty, rotx))

        # Direct calculation of angle and axis from A
        # Code adapted from thread response by Bruno Luong

        # Rotation axis u = [x, y, z]
        u = np.array([rotall[2, 1]-rotall[1, 2], rotall[0, 2]-rotall[2, 0], rotall[1, 0]-rotall[0, 1]])

        # Rotation sine and cosine
        c = np.trace(rotall) - 1
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


def timeseries_stats(s):
    """
    Compute useful stats for numpy timeseries vector
    Identify upper outliers in timeseries using UQ + 1.5 * IQR
    :param s: numpy vector of real values
    :return: tmean, thresh, nout, iout
        tmean = temporal mean of vector
        thresh = outlier threshold
        nout = number of outlier values
        iout = indices of outliers in timeseries
    """

    # Temporal mean of vector
    tmean = s.mean()

    # Find LQ and UQ of values
    lq, uq = np.percentile(s, (25, 75))

    # Outlier threshold = UQ + 1.5 (UQ - LQ) = 2.5 UQ - LQ
    thresh = 2.5 * uq - lq

    # Find indices of outliers in vector
    iout = np.array(np.where(s > thresh))

    # Count the outliers
    nout = iout.size

    # Percent of values that are outliers
    pout = nout / float(s.size) * 100.0

    return tmean, thresh, nout, pout, iout


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
