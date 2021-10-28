# Script to annotate intergenic insertions with promoter regions based on blastn of inferred TSS on Sa_MW2 (reference Prados, BMC Genomics 2016) using the internal reference as database

USAGE="Error! $0 takes one argument\\n
Usage: $0 <annotate intergenic directory>\\n
        $0 runs in a conda snippy environment"

if [ $# -eq 0 ]
then 
	echo -e $USAGE
	exit 1
fi

# Set variables
GROUP_DIR=$(readlink -e $1)
GROUP=$(basename $GROUP_DIR)
QUERY=$(readlink -e ~/saureus-general/ref-genomes/Sa_MW2/representative_promoter.fa)
BLAST=$(readlink -e ~/perl5/bin/run-blastn.sh)

# Initial statements
echo "This $0"
echo "This scripts annotate intergenic IS insertions with predicted promoters based on Prados, BMC Genomics"
echo "Annotating intergenic IS insertions in $GROUP"

# Check that the group directory contains intergenic mutations and exit if not
INTERGENIC=$(cat $GROUP_DIR/*/*.primary.intergenic.bed)
if [ -z $INTERGENIC ]
then
	echo "$GROUP doesn't have any intergenic insertions. Exiting"
	exit 1
fi

# Get promoters regions on the internal reference
cp -r $GROUP_DIR .

cd $GROUP

cd reference

mkdir promoters
cd promoters
cp ../ref.fa .
sh $BLAST $QUERY ref.fa promoters_ref_blastn
~/perl5/bin/blast2bed.sh promoters_ref_blastn.tab
mv promoters_ref_blastn.bed ref.promoters.bed

cd ../../

# Intersect with intergenic insertions
for i in $(cat isolates.txt)
do
	cd $i
	cut -f4-9 $i.primary.intergenic.bed > $i.intergenic.bed
	bedtools intersect -a ../reference/promoters/ref.promoters.bed -b $i.intergenic.bed -wa -wb > $i.primary.intergenic.promoters.bed
	rm $i.intergenic.bed
	cd ..
done

cd ..

# Final statements
echo "$0 done"
echo "You can find the annotated mutations in $GROUP/*/*/intergenic.promoters.bed"

