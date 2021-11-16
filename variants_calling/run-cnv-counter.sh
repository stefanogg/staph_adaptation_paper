#!/usr/bin/env bash

# Author: Stefano Giulieri
# Date: October 2021

# Shell script to run CNOGpro on the server
# Steps
# 1) Extract chromosome sequence from the reference (all strains within the group)
# 2) Get reads coverage as hits file
# 3) Run CNOGpro, normalise and run hmm. Save R object

# Run the script in conda snippy environment

# Functions
usage () {
	echo "USAGE"
	echo " $0 [options] <directory with bwa coverage output for the episode>"
	echo "OPTIONS:"
	echo " -f : force"
	echo " -R : resume"
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
RESUME=false
while getopts 'S:fR' option
do
	case $option in
		f) FORCE=true ;;
		R) RESUME=true ;;
	esac
done

# skip over the processed option to pass arguments
shift $((OPTIND-1))

# Create main directory
DIR=$(readlink -e $1)
GROUP_ID=$(basename $1)

# check if group directory exists and remove if true (later to be replaced by --force option)
if [ -d $GROUP_ID ]
then
	if $FORCE
	then
		echo "Directory $GROUP_ID exists. Removing"
		rm -r $GROUP_ID
	elif $RESUME
	then
		echo "Resuming analysis of $GROUP_ID"
	else
		echo "Directory $GROUP_ID exists. Use -f to remove directory. Exiting"
		exit 1
	fi
fi

# Create directory
mkdir $GROUP_ID
cd $GROUP_ID

# Copy reference and extract chromosome
mkdir reference
cd reference
cp $DIR/reference/*gbk ref.gbk
~/perl5/bin/split-gbk-chrom.py ref.gbk 
rm ref.gbk
any2fasta -u *.gbk > ref.fa
samtools faidx ref.fa
bedtools makewindows  -n 1 -g ref.fa.fai -n 1 > ref.bed
genbank2gff.pl *.gbk > ref.gff
cd ..

# Create internal reference folder
mkdir internal_ref
cp $DIR/internal_ref/internal_ref.txt internal_ref

# Extract list of isolates in the group
cut -f1 $DIR/*tab | sort -u > isolates.txt

# Run cnv-counter
for ISO in $(cat isolates.txt)
do
	~/perl5/bin/get-cnv-count.sh -S $SPLIT_DIR/$ISO $DIR/$ISO
done

# Copy internal ref files
INTERNAL_REF=$(cat internal_ref/internal_ref.txt)
cp $INTERNAL_REF/$INTERNAL_REF.cnv* internal_ref
if [ ! -f $INTERNAL_REF/$INTERNAL_REF.cnv.genes.bed ]
then
	touch internal_ref/$INTERNAL_REF.cnv.genes.bed
	touch internal_ref/$INTERNAL_REF.cnv.bed
fi

# Subtract cnv genes that are found in the internal reference
for ISO in $(cat isolates.txt)
do
	cd $ISO
	bedtools subtract -a $ISO.cnv.genes.bed -b ../internal_ref/$INTERNAL_REF.cnv.genes.bed > $ISO.mask.cnv.genes.bed
	bedtools subtract -a $ISO.cnv.bed -b ../internal_ref/$INTERNAL_REF.cnv.bed > $ISO.mask.cnv.bed
#	bedtools subtract -a $ISO.cnv.start.bed -b ../internal_ref/$INTERNAL_REF.cnv.start.bed > $ISO.mask.cnv.start.bed
#	bedtools subtract -a $ISO.cnv.stop.bed -b ../internal_ref/$INTERNAL_REF.cnv.stop.bed > $ISO.mask.cnv.stop.bed
	bedtools subtract -a $ISO.cnv.excision.pairs.bed -b ../internal_ref/$INTERNAL_REF.cnv.excision.pairs.bed > $ISO.mask.cnv.excision.pairs.bed
	cd ..
done

# Go back to main directory
cd ..
