# Shell pipeline to detect deletions based on regions of zero coverage on the closest available complete genome

echo "This is detect-deletions.sh"
echo "This script uses bedtools to detect deletions based on regions of zero coverage on the closest available complete genome"

# Usage

USAGE="$0 <grouped bwa directory>\\n
-f (force): removing existing directory"

if [ $# -eq 0 ]
then
	echo $USAGE
	exit 1
fi

# parse options
FORCE=false

while getopts 'f' option
do
	case $option in
		f) FORCE=true ;;
	esac
done

# skip to parse command line arguments
shift $((OPTIND-1))

# Parse command line arguments

BWA_DIR=$(readlink -e $1) # absolute path of bwa directory
PAIR_ID=$(basename $BWA_DIR) # pair id
LOG="detect-deletions.log"

echo "Copying $BWA_DIR bwa directory"

if [ -d $PAIR_ID ]
then
	if $FORCE
	then
		echo "Removing existing $PAIR_ID directory"
		rm -r $PAIR_ID
	else
		echo "Directory $PAIR_ID exists. Use -f to remove directory. Exiting"
               	exit 1
	fi
fi

# Create file structure

cp -r $BWA_DIR . # copy bwa directory
cd $PAIR_ID # move to copied bwa directory

echo "Running detect-deletions.sh on pair $PAIR_ID" 2>&1 | tee -a $LOG

# create internal reference folder
internal_ref=$(grep $PAIR_ID ../genetic_pairs_refgenomes.tab | cut -f3)
echo "Creating folder for internal reference $internal_ref coverage files" 2>&1 | tee -a $LOG

mkdir internal_ref
cd internal_ref
echo $internal_ref > internal_ref.txt
cd ..

# generate gff of the ref genome
cd reference
ref=$(basename *.gbk .gbk)
genbank2gff.pl $ref.gbk > $ref.gff
cd ..

# calculate coverage of all strains


for i in $(cut -f1 isolates.tab);
	do echo "Detecting deletions for strain $i"
	cd $i
	cd reference
	samtools faidx ref.fa
	bedtools makewindows -g ref.fa.fai -w 1000 > ref.bed
	cd ..
	echo "Running bedtools coverage with a window size of 1000 bp"
	bedtools coverage -a reference/ref.bed -b $i.bam > $i.coverage.bed
	awk '$7<1' $i.coverage.bed > $i.deletions.bed
	echo "Running bedtools genomecov with option -bga"
	bedtools genomecov -bga -ibam $i.bam | awk -v OFS='\t' '{print $0,$3-$2}' > $i.coverage.bedg
	awk '$4==0' $i.coverage.bedg > $i.deletions.bedg
	echo "Creating artemis userplot using samtools depth"
	samtools depth -aa $i.bam > $i.depth
	cut -f3 $i.depth > $i.userplot
	# move big coverage files to a separate folders (easier to remove later)
	echo "Moving large coverage files to a separate folder"
	mkdir coverage_files
	mv $i.coverage.bed* $i.depth coverage_files 
	cd ..
done 2>&1 | tee -a $LOG

# copy coverage files to internal reference folder
echo "Copying $internal_ref coverage files"
cd internal_ref
cp ../$internal_ref/$internal_ref.deletions.bed internal_ref.deletions.bed
cp ../$internal_ref/$internal_ref.deletions.bedg internal_ref.deletions.bedg
cd ..

# detect deletions that are not present on the internal reference
for i in $(cut -f1 isolates.tab);
	do echo "Detecting unique deletions in strain $i"
	cd $i
	bedtools subtract -a $i.deletions.bed -b ../internal_ref/internal_ref.deletions.bed > $i.mask.deletions.bed
	bedtools subtract -a $i.deletions.bedg -b ../internal_ref/internal_ref.deletions.bedg | awk -v OFS='\t' '{print $0,$3-$2}' > $i.mask.deletions.bedg
	cp ../reference/$ref.gff reference/ref.gff
	if [ -s $i.mask.deletions.bed ] 
	then
		del=$(wc -l $i.mask.deletions.bed | cut -f1 -d ' ')
		del_0=$(awk '$7==0' $i.mask.deletions.bed | wc -l)
		echo "Detected $del putative deleted 1000 bp-regions in strain $i, of which $del_0 have 0% coverage breadth" 
		bedtools merge -i $i.mask.deletions.bed -c 4,5,6,7 -o sum,sum,sum,mean > $i.merge.deletions.bed # alternative: bedtools cluster to keep original entries
		del_m=$(wc -l $i.merge.deletions.bed | cut -f1 -d ' ')
		del_m0=$(awk '$7==0' $i.merge.deletions.bed | wc -l)
		echo "After merging detected $del_m putative deleted regions in strain $i, of which $del_m0 have 0% coverage breadth"
	else
		echo "No uniquely deleted regions detected in strain $i (bed file). Is $i the internal reference?"
	fi
	if [ -s $i.mask.deletions.bedg ]
	then
		del_g=$(wc -l $i.mask.deletions.bedg | cut -f1 -d ' ')
		del_gg=$(awk '$6>500' $i.mask.deletions.bedg | wc -l)
		echo "Detected $del_g deleted regions in strain $i, of which $del_gg are > 500 bp"
		bedtools intersect -a $i.mask.deletions.bedg -b reference/ref.gff -wa -wb > $i.mask.deletions.gff
	else
		echo "No uniquely deletions detected in strain $i (bedg file). Is $i the internal reference?"
	fi
	# Move deletions file to a separate folder
	echo "Moving deletions files to a separate folder. Only unique deletions files (*.mask.deletions.*) are left in the main directory of strain $i"
	mkdir deletions_files
	mv $i.deletions.bed* deletions_files
	cd ..
done 2>&1 | tee -a $LOG

cd .. 
