#!/bin/bash

# Download FASTQ files from core's FTP site using sftp. Must be run
# interactively because a password is required. Downloads to a file called
# FastQ.
#
# Usage:
#
# bash fastq-download.sh remote local
#
# For example:
#
# bash fastq-download.sh NGS-2017/171020_700819F_0579_ACAPCDACXX-YG-PYT-Flowcell-A/FastQ /project2/gilad/singlecell-qtl/fastq

remote=$1
local=$2

echo -e "Remote:\t$remote"
echo -e "Local:\t$local"

mkdir -p $local

sftp -r gilad@fgfftp.uchicago.edu:/Genomics/$remote $local
