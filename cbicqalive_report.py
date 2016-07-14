#!/usr/bin/env python3
#
# Create daily QA HTML report
#
# USAGE : cbicqa_report.py <QA Directory>
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 09/25/2013 JMT From scratch
#          10/23/2013 JMT Add com external call
#          10/24/2013 JMT Move stats calcs to new cbicqa_stats.py
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
import string
import argparse
import numpy as np
from pylab import *

# Define template
TEMPLATE_FORMAT = """
<html>

<head>
<STYLE TYPE="text/css">
BODY {
  font-family    : arial, sans-serif;
}
td, th {
  padding-left   : 10px;
  padding-right  : 10px;
  padding-top    : 0px;
  padding-bottom : 0px;
  text-align     : "left";
}
</STYLE>
</head>

<body>

<h1 style="background-color:#E0E0FF">CBIC In Vivo Quality Assurance</h1>

<h3> Dataset Information </h3>
<table>
    <tr>
        <td> QA Directory </td>
        <td> $qa_dir_abs </td>
    </tr>
    <tr>
        <td> Repetition Time (seconds) </td>
        <td> $TR_secs </td>
    </tr>
    <tr>
        <td> Image Volumes </td>
        <td> $N_vols </td>
    </tr>

</table>

<h3> Summary Statistics </h3>
<table>
    <tr>
        <td> <b> Parameter </b> </td>
        <td> <b> Mean </b> </td>
        <td> <b> Threshold </b> </td>
        <td> <b> Outliers </b> </td>
        <td> </td>
    </tr>
    <tr>
        <td> sMean Brain Signal </td>
        <td> $brain_tmean </td>
        <td> $brain_thresh </td>
        <td> $brain_nout </td>
        <td> $brain_pout% </td>
    </tr>
    <tr>
        <td> sMean Nyquist Ghost Signal</td>
        <td> $ghost_tmean </td>
        <td> $ghost_thresh </td>
        <td> $ghost_nout </td>
        <td> $ghost_pout% </td>
    </tr>
    <tr>
        <td> sMean Air Signal</td>
        <td> $air_tmean </td>
        <td> $air_thresh </td>
        <td> $air_nout </td>
        <td> $air_pout% </td>
    </tr>
    <tr>
        <td> DVARS</td>
        <td> $dvars_tmean </td>
        <td> $dvars_thresh </td>
        <td> $dvars_nout </td>
        <td> $dvars_pout% </td>
    </tr>
    <tr>
        <td> F-F Displacement (microns)</td>
        <td> $dd_um_tmean </td>
        <td> $dd_um_thresh </td>
        <td> $dd_um_nout </td>
        <td> $dd_um_pout% </td>
    </tr>
    <tr>
        <td> F-F Rotation (mdeg)</td>
        <td> $dr_mdeg_tmean </td>
        <td> $dr_mdeg_thresh </td>
        <td> $dr_mdeg_nout </td>
        <td> $dr_mdeg_pout% </td>
    </tr>
</table>

<h3> Temporal Summary Images </h3>
<table>
    <tr>
        <td> <b>SNFR</b><br> <img src=qa_tsnfr_ortho.png /> </td>
    </tr>
    <tr>
        <td> <br><b>Mean Signal</b><br> <img src=qa_tmean_ortho.png /> </td>
    </tr>
    </tr>
        <td> <br><b>Fluctuation Noise SD</b><br> <img src=qa_tsd_ortho.png /> </td>
    </tr>
    <tr>
        <td> <br><b>Region Mask</b><br> <img src=qa_mask_ortho.png /> </td>
    </tr>
</table>

<h3> Signal Timeseries </h3>
<table>
    <tr>
        <td> <img src=qa_roi_timeseries.png /> </td>
    <tr>
</table>

<h3> Motion Timeseries </h3>
<table>
    <tr>
        <td> <img src=qa_motion_timeseries.png /><br> </td>
    </tr>
</table>

</body>
"""


# Main function
def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='QA reporting for in vivo fMRI timeseries')
    parser.add_argument('-i', '--qa_dir', help="CBICQAlive directory (*.qalive)")

    # Parse command line arguments
    args = parser.parse_args()
    qa_dir = args.qa_dir
    
    print('  Creating in vivo QA report for ' + qa_dir)

    # Determine full path for QA directory
    qa_dir_abs = os.path.abspath(qa_dir)

    # Load dataset info from QA directory
    info_fname = os.path.join(qa_dir, 'qa_info.txt')
    if not os.path.isfile(info_fname):
        print(info_fname + ' does not exist - exiting')
        sys.exit(0)

    info = np.genfromtxt(info_fname, delimiter=',')
    TR_secs = info[0,1]
    N_vols = info[1,1]

    # Load summary statistics from CSV file in QA directory
    stats_csv_fname = os.path.join(qa_dir, 'qa_stats.csv')
    if not os.path.isfile(stats_csv_fname):
        print(stats_csv_fname + ' does not exist - exiting')
        sys.exit(0)

    stats = np.genfromtxt(stats_csv_fname, delimiter=',', usecols=(1,2,3,4))

    #
    # HTML report generation
    #

    # Create substitution dictionary for HTML report
    qa_dict = dict([
        ('qa_dir_abs',     "%s"    % qa_dir_abs),
        ('TR_secs',        "%0.3f" % TR_secs),
        ('N_vols',         "%d"    % N_vols),
        ('brain_tmean',    "%0.1f" % stats[0, 0]),
        ('brain_thresh',   "%0.1f" % stats[0, 1]),
        ('brain_nout',     "%d"    % stats[0, 2]),
        ('brain_pout',     "%0.1f" % stats[0, 3]),
        ('ghost_tmean',    "%0.1f" % stats[1, 0]),
        ('ghost_thresh',   "%0.1f" % stats[1, 1]),
        ('ghost_nout',     "%d"    % stats[1, 2]),
        ('ghost_pout',     "%0.1f" % stats[1, 3]),
        ('air_tmean',      "%0.1f" % stats[2, 0]),
        ('air_thresh',     "%0.1f" % stats[2, 1]),
        ('air_nout',       "%d"    % stats[2, 2]),
        ('air_pout',       "%0.1f" % stats[2, 3]),
        ('dvars_tmean',    "%0.1f" % stats[3, 0]),
        ('dvars_thresh',   "%0.1f" % stats[3, 1]),
        ('dvars_nout',     "%d"    % stats[3, 2]),
        ('dvars_pout',     "%0.1f" % stats[3, 3]),
        ('dd_um_tmean',    "%0.1f" % stats[4, 0]),
        ('dd_um_thresh',   "%0.1f" % stats[4, 1]),
        ('dd_um_nout',     "%d"    % stats[4, 2]),
        ('dd_um_pout',     "%0.1f" % stats[4, 3]),
        ('dr_mdeg_tmean',  "%0.1f" % stats[5, 0]),
        ('dr_mdeg_thresh', "%0.1f" % stats[5, 1]),
        ('dr_mdeg_nout',   "%d"    % stats[5, 2]),
        ('dr_mdeg_pout',   "%0.1f" % stats[5, 3]),
    ])

    # Generate HTML report from template (see above)
    TEMPLATE = string.Template(TEMPLATE_FORMAT)
    html_data = TEMPLATE.safe_substitute(qa_dict)
    
    # Write HTML report page
    qa_report_file = os.path.join(qa_dir, 'index.html')
    open(qa_report_file, "w").write(html_data)

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
