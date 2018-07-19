#!/usr/bin/env python
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
# The MIT License (MIT)
#
# Copyright (c) 2018 Mike Tyszka
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
            <td> <b> Median tSFNR </b>
            <td> <b> $tSFNR </b>
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
            <td> <b> Percent Outliers 
        </tr>
            <td> Signal  
            <td> $signal_tmean 
            <td> $signal_thresh 
            <td> $signal_pout%  
        </tr>
        <tr>
            <td> Nyquist Ghost  
            <td> $ghost_tmean 
            <td> $ghost_thresh 
            <td> $ghost_pout% 
        </tr>
        <tr>
            <td> Air 
            <td> $air_tmean 
            <td> $air_thresh 
            <td> $air_pout% 
        </tr>
        <tr>
            <td> DVARS
            <td> $dvars_tmean 
            <td> $dvars_thresh 
            <td> $dvars_pout% 
        </tr>
        <tr>
            <td> F-F Displacement (microns)
            <td> $dd_um_tmean 
            <td> $dd_um_thresh 
            <td> $dd_um_pout% 
        </tr>
        <tr>
            <td> F-F Rotation (mdeg)
            <td> $dr_mdeg_tmean 
            <td> $dr_mdeg_thresh 
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
            <td> <br><b>Region Mask</b><br> <img src=qc_labels_ortho.png /> 
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

    # Load dataset info and stats from QC directory
    info_fname = os.path.join(qc_dir, 'qc_info.csv')
    info = genfromtxt(info_fname, delimiter=',')

    # Create substitution dictionary for HTML report
    qc_dict = dict([
        ('qc_dir_abs',     "%s"    % qc_dir_abs),
        ('TR_secs',        "%0.3f" % info[1]),
        ('N_vols',         "%d"    % info[2]),
        ('tSFNR',          "%0.1f" % info[3]),
        ('signal_tmean',    "%0.1f" % info[4]),
        ('signal_thresh',   "%0.1f" % info[5]),
        ('signal_pout',     "%0.1f" % info[6]),
        ('ghost_tmean',    "%0.1f" % info[7]),
        ('ghost_thresh',   "%0.1f" % info[8]),
        ('ghost_pout',     "%0.1f" % info[9]),
        ('air_tmean',      "%0.1f" % info[10]),
        ('air_thresh',     "%0.1f" % info[11]),
        ('air_pout',       "%0.1f" % info[12]),
        ('dvars_tmean',    "%0.1f" % info[13]),
        ('dvars_thresh',   "%0.1f" % info[14]),
        ('dvars_pout',     "%0.1f" % info[15]),
        ('dd_um_tmean',    "%0.1f" % info[16]),
        ('dd_um_thresh',   "%0.1f" % info[17]),
        ('dd_um_pout',     "%0.1f" % info[18]),
        ('dr_mdeg_tmean',  "%0.1f" % info[19]),
        ('dr_mdeg_thresh', "%0.1f" % info[20]),
        ('dr_mdeg_pout',   "%0.1f" % info[21]),
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
