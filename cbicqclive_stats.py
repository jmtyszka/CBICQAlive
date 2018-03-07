#!/usr/bin/env python
#
# Derive various descriptive stats from the preprocessed QA data
# - detrended time series parameters
# - spike counts
# - center of mass
# - apparent motion
#
# USAGE : cbicqc_stats.py <QA Directory>
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 04/01/2013 JMT From scratch
#          10/24/2013 JMT Expand to calculate all stats
#          03/03/2014 JMT Add scanner name argument
#          2017-07-31 JMT Refactor quality assurance to quality control
#
# This file is part of CBICQC.
#
#    CBICQC is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    CBICQC is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#   along with CBICQC.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2013-2017 California Institute of Technology.

import os
import sys
import argparse
import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# Main function
def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='CBICQClive statistics')
    parser.add_argument('-i', '--qc_dir', help="CBICQClive directory (*.qalive)")

    # Parse command line arguments
    args = parser.parse_args()
    qc_dir = args.qc_dir

    print('  Calculating statistics for ' + qc_dir)
    
    # Load ROI timeseries
    print('    Loading ROI timeseries')
    qc_ts_fname = os.path.join(qc_dir, 'qc_roi_ts.txt')
    if not os.path.isfile(qc_ts_fname):
        print(qc_ts_fname + ' does not exist - exiting')
        sys.exit(0)

    x = np.loadtxt(qc_ts_fname)

    # Parse timeseries into vectors
    # Three timeseries: Signal, Nyquist, Air
    signal_t = x[:, 0]
    ghost_t = x[:, 1]
    air_t   = x[:, 2]

    # Timeseries vector
    nv = len(signal_t)
    v = np.linspace(0, nv-1, nv)

    # Load DVARS timeseries
    print('    Loading DVARS timeseries')
    qc_dvars_fname = os.path.join(qc_dir, 'qc_dvars_ts.txt')
    if not os.path.isfile(qc_dvars_fname):
        print(qc_dvars_fname + ' does not exist - exiting')
        sys.exit(0)

    dvars_t = np.loadtxt(qc_dvars_fname)

    #
    # Analyze motion parameters
    #
    print('    Analyzing motion parameters')
    qc_mcf_parfile = os.path.join(qc_dir, 'qc_mcf.par')

    if not os.path.isfile(qc_mcf_parfile):
        print(qc_mcf_parfile + ' does not exist - exiting')
        sys.exit(0)

    # Read file into 2D array
    try:
        mcf_pars = np.loadtxt(qc_mcf_parfile)
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
    signal_tmean, signal_thresh, signal_nout, signal_pout, _ = timeseries_stats(signal_t)
    ghost_tmean, ghost_thresh, ghost_nout, ghost_pout, _ = timeseries_stats(ghost_t)
    air_tmean, air_thresh, air_nout, air_pout, _ = timeseries_stats(air_t)
    dvars_tmean, dvars_thresh, dvars_nout, dvars_pout, _ = timeseries_stats(dvars_t)
    dd_um_tmean, dd_um_thresh, dd_um_nout, dd_um_pout, _ = timeseries_stats(dd_um)
    dr_mdeg_tmean, dr_mdeg_thresh, dr_mdeg_nout, dr_mdeg_pout, _ = timeseries_stats(dr_mdeg)

    #
    # Output summary statistics to CSV file
    #
    stats_csv_fname = os.path.join(qc_dir, 'qc_stats.csv')
    fd = open(stats_csv_fname, "w")

    fd.write("Signal, %0.1f, %0.1f, %d, %0.1f\n" % (signal_tmean, signal_thresh, signal_nout, signal_pout))
    fd.write("Ghost, %0.1f, %0.1f, %d, %0.1f\n" % (ghost_tmean, ghost_thresh, ghost_nout, ghost_pout))
    fd.write("Air, %0.1f, %0.1f, %d, %0.1f\n" % (air_tmean, air_thresh, air_nout, air_pout))
    fd.write("DVARS, %0.1f, %0.1f, %d, %0.1f\n" % (dvars_tmean, dvars_thresh, dvars_nout, dvars_pout))
    fd.write("FF_Disp_um, %0.1f, %0.1f, %d, %0.1f\n" % (dd_um_tmean, dd_um_thresh, dd_um_nout, dd_um_pout))
    fd.write("FF_Rot_mdeg, %0.1f, %0.1f, %d, %0.1f\n" % (dr_mdeg_tmean, dr_mdeg_thresh, dr_mdeg_nout, dr_mdeg_pout))

    fd.close()

    # ROI timeseries figure

    fig = plt.figure(figsize=(10, 10))

    plt.subplot(411)
    plt.plot(v, signal_t)
    plt.axhline(y=signal_thresh, xmin=0, xmax=1, linestyle=':', linewidth=2, color='r')
    plt.title("sMean Signal", x=0.5, y=0.8)

    plt.subplot(412)
    plt.plot(v, ghost_t)
    plt.axhline(y=ghost_thresh, xmin=0, xmax=1, linestyle=':', linewidth=2, color='r')
    plt.title("sMean Ghost", x=0.5, y=0.8)

    plt.subplot(413)
    plt.plot(v, air_t)
    plt.axhline(y=air_thresh, xmin=0, xmax=1, linestyle=':', linewidth=2, color='r')
    plt.title("sMean Air", x=0.5, y=0.8)

    plt.subplot(414)
    plt.plot(v, dvars_t)
    plt.axhline(y=dvars_thresh, xmin=0, xmax=1, linestyle=':', linewidth=2, color='r')
    plt.title("DVARS", x=0.5, y=0.8)

    plt.savefig(os.path.join(qc_dir, 'qc_roi_timeseries.png'), dpi=72, bbox_inches='tight')

    # Motion timeseries figure

    fig.clf()

    plt.subplot(411)
    plt.plot(v, x_um, 'r', v, y_um, 'g', v, z_um, 'b')
    plt.title("Displacement (microns)", x=0.5, y=0.8)

    plt.subplot(412)
    plt.plot(v, rx_mdeg, 'r', v, ry_mdeg, 'g', v, rz_mdeg, 'b')
    plt.title("Rotation (mdeg)", x=0.5, y=0.8)

    plt.subplot(413)
    plt.plot(v, dd_um)
    plt.axhline(y=dd_um_thresh, xmin=0, xmax=1, linestyle=':', linewidth=2, color='r')
    plt.title("Total F-F Displacement (microns)", x=0.5, y=0.8)

    plt.subplot(414)
    plt.plot(v, dr_mdeg)
    plt.axhline(y=dr_mdeg_thresh, xmin=0, xmax=1, linestyle=':', linewidth=2, color='r')
    plt.title("Total F-F Rotation (mdeg)", x=0.5, y=0.8)

    plt.savefig(os.path.join(qc_dir, 'qc_motion_timeseries.png'), dpi=72, bbox_inches='tight')

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

    :param s: numpy float array
        Timeseries

    :return:
        tmean: float
            Timeseries mean
        thresh: float
            Outlier threshold
        nout: int
            Number of outlier values
        iout: numpy int array
            Indices of outliers in timeseries

    """

    # Mean of timeseries
    tmean = s.mean()

    # Find LQ and UQ of timeseries
    lq, uq = np.percentile(s, (25, 75))

    # Outlier threshold = UQ + 1.5 (UQ - LQ)
    thresh = uq + 1.5 * (uq - lq)

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
