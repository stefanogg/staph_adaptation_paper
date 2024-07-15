#!/bin/bash
#
# This is a bash script to run provean given: 1) the locus tag of the gene; 2) the list of non-synonymous substitutions
# It also requires a multifasta files with the amino-acid sequences and a tab-separated file with the genes and substitutions
# It needs provean and seqkit: run in a conda base environment"
#
# Author: Stefano Giulieri
# Date: April 2022

# Functions
# Usage
usage () {
        echo "USAGE"
        echo " $0 [options] <LOCUS-TAG>"
        echo "OPTIONS:"
        echo " -A : fasta file with amino-acid sequences (default: mutated_genes.faa)"
        echo " -I : input file (default: provean_in.tab)"
        echo " -f : force (default: false)"
        echo "Please run in a conda base environment"
}

# Print usage and exit if no arguments
if [ $# -eq 0 ]
then
        usage
        exit 1
fi

# PROVEAN location
PROVEAN="/home/giulieris/perl5/bin/provean/provean.sh"

# Parse options
FASTA="mutated_genes.faa"
IN="provean_in.tab"
FORCE=false

while getopts 'A:I:f' option
do
        case $option in
                A) FASTA=$OPTARG ;;
                I) IN=$OPTARG ;;
                f) FORCE=true ;;
        esac
done

# Get absolute path to fhe files
FASTA=$(readlink -e $FASTA)
IN=$(readlink -e $IN)

# skip to parse command line arguments
shift $((OPTIND-1))

# Get positional arguments
LOCUS_TAG=$1

# Starting statements
echo "Running $0"
echo "The fasta file with the amino-acid sequences is $FASTA"
echo "The input file is $IN"

# Parse input file and prepare file structure
echo "Extracting variants for $LOCUS_TAG"

if [ -d $LOCUS_TAG ]
then
        if $FORCE
        then
                echo "Directory $LOCUS_TAG exists. Removing"
                rm -r $LOCUS_TAG
        else
                echo "Directory $LOCUS_TAG exists. Use -f to remove directory. Exiting"
                exit 1
        fi
fi

mkdir $LOCUS_TAG
cd $LOCUS_TAG

seqkit grep -p $LOCUS_TAG $FASTA > $LOCUS_TAG.faa
grep $LOCUS_TAG $IN | cut -f2 > $LOCUS_TAG.var

# Run provean
echo "Running PROVEAN"
$PROVEAN -q $LOCUS_TAG.faa -v $LOCUS_TAG.var --save_supporting_set $LOCUS_TAG.sss > $LOCUS_TAG.out

# Final statement
cd ..
echo "$0 done"
