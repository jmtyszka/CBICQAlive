# CBICQAlive
MRI QA script for in vivo EPI data

## Dependencies:
* Python 3.4
* FSL 5x

## Installation

Clone the repository into a suitable location (eg /usr/local) and add this directory to your BASH path

## Usage

cbicqalive \<summary statistics file prefix\> \<list of 4D EPI compressed Nifti images\>

## Example

If you have a dataset consisting of a group of patients with multiple 4D EPI fMRI timeseries two levels down in the directory tree, then the script can produce HTML in vivo QA reports for each 4D EPI image as follows:

cbicqalive summary_stats \*\/\*\/epi.nii.gz

A summary table of dataset information and stats results will be created in summary_stats.csv (the first command line argument)
epi.qalive folders containing intermediate files and the HTML report (index.html) will be created at the same level as the epi.nii.gz files.
