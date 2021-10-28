# This is a shell script to get intergenic regions from a list of IS insertions, annoted them from the gff of the reference and get the fasta sequences
# It takes one positional argument: the directory with the run-annotate-splitters.sh output

# Usage

USAGE="$0 <annotate-splitters directory>"

if [ $# -eq 0 ]
then
	echo $USAGE
	exit 1
fi


ANNOT_DIR=$(readlink -e $1) # grouped directory with output from run-annotate-splitters.sh
PAIR_ID=$(basename $ANNOT_DIR) # pair id


# create new directory
[[ -d $PAIR_ID ]] && (echo "Removing existing directory $PAIR_ID"; rm -r $PAIR_ID)
mkdir $PAIR_ID
cd $PAIR_ID

# Get intergenic regions of the reference

cp -r $ANNOT_DIR/reference .

cd reference
grep CDS ref.gff > ref.cds.gff
bedtools makewindows -g ref.fa.fai -n 1 | bedtools subtract -a stdin -b ref.cds.gff > ref.intergenic.bed

cd ..

# Get interegenic regions of insertions

for i in $(cat $ANNOT_DIR/isolates.txt)
do
	# check if isolates has splitters	
	bed="$ANNOT_DIR/$i/$i.primary.sa.bed"
	if [ ! -f $bed ]
	then
		echo "Strain $i doesn't have new IS insertions. Interrupting the annotation"	
		continue
	fi

	echo $i >> isolates.txt
	mkdir $i
	cd $i
	grep primary $bed > $i.primary.bed
	bedtools intersect -a ../reference/ref.intergenic.bed -b $i.primary.bed -wa -wb > $i.primary.intergenic.bed
	bedtools getfasta -fi ../reference/ref.fa -bed $i.primary.intergenic.bed > $i.primary.intergenic.fa
	bedtools window -a ../reference/ref.cds.gff -b $i.primary.intergenic.bed -w 1 >  $i.primary.intergenic.gff
	cd ..
done


