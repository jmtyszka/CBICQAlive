#!/usr/bin/env python3
#
# Create daily QC HTML report
#
# USAGE : cbicqc_report.py <QA Directory>
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 09/25/2013 JMT From scratch
#          10/23/2013 JMT Add com external call
#          10/24/2013 JMT Move stats calcs to new cbicqc_stats.py
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
# Copyright 2013-2014 California Institute of Technology.

import os
import string
import argparse
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

<h1 style="background-color:#E0E0FF">CBIC In Vivo Quality Control</h1>

<div>
    <table>
        <tr>
            <td> QC Directory 
            <td> $qc_dir_abs 
        </tr>
        <tr>
            <td> Repetition Time (seconds) 
            <td> $TR_secs 
        </tr>
        <tr>
            <td> Image Volumes 
            <td> $N_vols 
        </tr>
        <tr style="background-color:#AFFF9F">
            <td> <b> Brain tSFNR </b>
            <td> <b> $tSFNR_brain </b>
        </tr>    
    </table>
</div>

<br>

<div>
    <table>
        <tr>
            <td> <b> Parameter </b> 
            <td> <b> tMean </b> 
            <td> <b> Threshold </b> 
            <td> <b> Outliers</b> 
            <td> <b> Outlier Fraction 
        </tr>
            <td> Brain  
            <td> $brain_tmean 
            <td> $brain_thresh 
            <td> $brain_nout 
            <td> $brain_pout%  
        </tr>
        <tr>
            <td> Nyquist Ghost  
            <td> $ghost_tmean 
            <td> $ghost_thresh 
            <td> $ghost_nout 
            <td> $ghost_pout% 
        </tr>
        <tr>
            <td> Air 
            <td> $air_tmean 
            <td> $air_thresh 
            <td> $air_nout 
            <td> $air_pout% 
        </tr>
        <tr>
            <td> DVARS
            <td> $dvars_tmean 
            <td> $dvars_thresh 
            <td> $dvars_nout 
            <td> $dvars_pout% 
        </tr>
        <tr>
            <td> F-F Displacement (microns)
            <td> $dd_um_tmean 
            <td> $dd_um_thresh 
            <td> $dd_um_nout 
            <td> $dd_um_pout% 
        </tr>
        <tr>
            <td> F-F Rotation (mdeg)
            <td> $dr_mdeg_tmean 
            <td> $dr_mdeg_thresh 
            <td> $dr_mdeg_nout 
            <td> $dr_mdeg_pout% 
        </tr>
    </table>
</div>

<br>

<div>
    <table>
        <tr>
            <td> <br><b>Motion Timeseries</b><br> <img src=qc_motion_timeseries.png /> 
            <td> <br><b>ROI Timeseries</b><br> <img src=qc_roi_timeseries.png /> 
        <tr>
    </table>
</div>

<div>
    <table>
        <tr>
            <td> <br><b>Temporal Mean Signal</b><br> <img src=qc_tmean_ortho.png /> 
            <td> <br><b>Fluctuation Noise SD</b><br> <img src=qc_tsd_ortho.png /> 
        </tr>
        <tr>
            <td> <b>Temporal Signal-to-Fluctuation-Noise Ratio (SFNR)</b><br> <img src=qc_tsfnr_ortho.png /> 
            <td> <br><b>Region Mask</b><br> <img src=qc_mask_ortho.png /> 
        </tr>
    </table>
</div>

</body>
"""


# Main function
def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='QC reporting for in vivo fMRI timeseries')
    parser.add_argument('-i', '--qc_dir', help="CBICQClive directory (*.qclive)")

    # Parse command line arguments
    args = parser.parse_args()
    qc_dir = args.qc_dir
    
    print('  Creating in vivo QC report for ' + qc_dir)

    # Determine full path for QC directory
    qc_dir_abs = os.path.abspath(qc_dir)

    # Load dataset info from QC directory
    info_fname = os.path.join(qc_dir, 'qc_info.txt')
    if not os.path.isfile(info_fname):
        print(info_fname + ' does not exist - exiting')
        sys.exit(0)

    info = np.genfromtxt(info_fname, delimiter=',')
    TR_secs = info[0,1]
    N_vols = info[1,1]

    # Load summary statistics from CSV file in QC directory
    stats_csv_fname = os.path.join(qc_dir, 'qc_stats.csv')
    if not os.path.isfile(stats_csv_fname):
        print(stats_csv_fname + ' does not exist - exiting')
        sys.exit(0)

    stats = np.genfromtxt(stats_csv_fname, delimiter=',', usecols=(1,2,3,4))

    # Load sMean tSFNR results from file
    sfnr_fname = os.path.join(qc_dir, 'qc_tsfnr_brain.txt')
    if not os.path.isfile(sfnr_fname):
        print(sfnr_fname + ' does not exist - exiting')
        sys.exit(0)
    tSFNR_brain = np.genfromtxt(sfnr_fname)

    #
    # HTML report generation
    #

    # Create substitution dictionary for HTML report
    qc_dict = dict([
        ('qc_dir_abs',     "%s"    % qc_dir_abs),
        ('TR_secs',        "%0.3f" % TR_secs),
        ('N_vols',         "%d"    % N_vols),
        ('tSFNR_brain',    "%0.1f" % tSFNR_brain),
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
    html_data = TEMPLATE.safe_substitute(qc_dict)
    
    # Write HTML report page
    qc_report_html = os.path.join(qc_dir, 'index.html')
    open(qc_report_html, "w").write(html_data)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
