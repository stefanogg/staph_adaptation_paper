#!/usr/bin/env bash

# Author: Stefano Giulieri
# Date: October 2021

# Shell script to run CNOGpro on the server
# Steps
# 1) Create isolate directory
# 2) Copy reference
# 3) Get reads coverage as hits file
# 4) Run CNOGpro, normalise and run hmm. Save R object

# Run the script in conda snippy environment

# Functions
usage () {
	echo "USAGE"
	echo " $0 [options] <directory with bwa coverage output for the isolate>"
	echo "OPTIONS:"
	echo " -S : directory with split break output for the isolate"
	echo " -f : force"
	echo "Please run this script in conda snippy environment"
}

# Print usage and exit if no arguments
if [ $# -eq 0 ]
then
    	usage
	exit 1
fi

# Parse options
FORCE=false
while getopts 'S:f' option
do
	case $option in
		S) SPLIT_DIR=$OPTARG ;;
		f) FORCE=true ;;
	esac
done

# skip over the processed option to pass arguments
shift $((OPTIND-1))

# Create main directory
DIR=$(readlink -e $1)
ISO=$(basename $DIR)
SPLIT_DIR=$(readlink -e $SPLIT_DIR)

# Opening statement
echo "This $0"
echo "Getting copy number variants for $ISO"

# check if group directory exists and remove if true (later to be replaced by --force option)
if [ -d $ISO ]
then
	if $FORCE
	then
		echo "Directory $ISO exists. Removing"
		rm -r $ISO
	else
		echo "Directory $ISO exists. Use -f to remove directory. Exiting"
		exit 1
	fi
fi

# Create directory
mkdir $ISO
cd $ISO

# Copy reference and extract chromosome
mkdir reference
cd reference
cp ../../reference/*gbk ref.gbk
cp ../../reference/ref.* .
cd ..

# Generate hits file
samtools view -L reference/ref.bed $DIR/$ISO.bam | cut -f3,4 > $ISO.hits

# Run CNOGpro
eval "$(conda shell.bash hook)"
conda activate R
Rscript --default-packages=CNOGpro,utils,stats ~/R_functions/run-CNOGpro.R $ISO 2>&1 | tee -a  $ISO.log
# Get overlap with split reads
if [ -z $SPLIT_DIR ]
then
	echo "No split-break output available for $ISO. Exiting"
	cd ..
	exit 1
fi

eval "$(conda shell.bash hook)"
conda activate snippy
bedtools merge -i $ISO.cnv.genes.bed -d 1 | bedtools window -l 10000 -r 0 -a stdin -b $SPLIT_DIR/$ISO.start.bed > $ISO.cnv.start.bed
# optional keep reads on the left the duplication only
# awk '$10<$2'
bedtools merge -i $ISO.cnv.genes.bed -d 1 | bedtools window -l 0 -r 10000 -a stdin -b $SPLIT_DIR/$ISO.stop.bed > $ISO.cnv.stop.bed
bedtools intersect -wa -wb -a $ISO.cnv.start.bed -b $ISO.cnv.stop.bed > $ISO.cnv.excision.pairs.bed

# Go back to main directory
cd ..
