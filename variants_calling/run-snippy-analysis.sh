#!/usr/bin/env bash
# Pipeline to run snippy on patient-episode isolates

# Author: Stefano Giulieri
# Date: July 2020

USAGE="Usage: run-snippy-analysis.sh <patient-episode-id>"

# ---- Check arguments ---------------------------

if [ $# == 0 ]; then
    echo -e $USAGE
    exit 1
fi

# Set variables
peid=$1

# ----- Run script --------------------------------------------
echo "Running snippy on	patient-episode	$peid"

cd $peid

# Get reference gbk
id=$(cut -f1 reference.tab)
dir=$(cut -f3 reference.tab)
gbk="$dir/$id.gbk"
# echo $gbk
echo "Reference is $gbk"

mkdir reference
cp $gbk reference


# Run snippy on pe reads
echo "Running snippy on pair-end reads"
cat reads.tab | while read f1 f2 f3;
	do snippy --outdir $f1 --ref reference/*gbk --R1 $f2 --R2 $f3
done

# Run snippy on shredded assembly
echo "Running snippy on the shredded assembly"
cut -f1,3 assembly.tab | while read f1 f2;
	do a="$f2/$f1.fna"
	snippy --outdir $f1 --ref reference/*gbk --ctgs $a
done # need to remove f2 from assembly.tab because read doesn't skip empty columns (there is probably a way around this - to be enquired later)


# Count number of variants
for i in $(cut -f1 reads.tab ) $(cut -f1 assembly.tab ); do (echo $i; grep VariantTotal $i/snps.txt | cut -f2) | paste - -; done > all_variants.count.tab


# Run snippy-core
echo "Running snippy-core"
snippy-core --ref reference/*gbk $(cut -f1 reads.tab) $(cut -f1	assembly.tab)

# Run snippy-coreness
echo "Running snippy-coreness"
snippy-coreness core.full.aln > coreness.tab

cd ..
