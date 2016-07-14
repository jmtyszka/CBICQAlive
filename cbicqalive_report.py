#!/opt/local/bin/python
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

import sys
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
  font-family    : sans-serif;
}
td {
  padding-left   : 10px;
  padding-right  : 10px;
  padding-top    : 0px;
  padding-bottom : 0px;
  vertical-align : top;
}
</STYLE>
</head>

<body>

<h1 style="background-color:#E0E0FF">CBIC Daily QA Report</h1>

<table>
<tr>

<!-- Scanner and acquisition info -->
<td>
  <table>
  <tr> <td> <b>Dummy</b> <td bgcolor="#E0FFE0"> <b>XXXXXX</b> </tr>
  </table>
</td>

<br><br>

<table>

<tr>
<td> <h3>Temporal Summary Images and Masks</h3>
<td> <h3>Signal, Drift and Noise</h3>
</tr>

<tr>

<!-- Temporal summary images and masks -->
<td valign="top">
<b>tSNR</b><br> <img src=qa_tsnr_ortho.png /><br><br>
<b>tMean</b><br> <img src=qa_tmean_ortho.png /><br><br>
<b>tSD</b><br> <img src=qa_tsd_ortho.png /><br><br>
<b>Region Mask</b><br> <img src=qa_mask_ortho.png /><br><br>

<!-- Motion timeseries -->
<td valign="top">
<img src=qa_mcf_trans.png /><br><br>
<img src=qa_mcf_rot.png /><br><br>

</tr>

</table>

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

    # Load QA acquisition info
    try:
		    qa_info_file = os.path.join(qa_dir, 'qa_info.txt')
        # Load metadata from info file
    except:
        print('*** No metadata detected in qa_info.txt - skipping')

    #
    # HTML report generation
    #

    # Create substitution dictionary for HTML report
    qa_dict = dict([
      ('phantom_mean',   "%0.1f" % (123.456)),
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
