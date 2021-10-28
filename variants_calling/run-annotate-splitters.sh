# Shell wrapper for the annotate-splitters script
# It has a single argument: the directory with the split break output for the patient episode or pair
# Note that it needs to be run in the conda snippy environment

# Usage
usage () {
	echo "USAGE"
        echo "  $0 [options] <split-break directory>"
        echo "OPTIONS"
        echo "  -f (force): remove existing directory"
}

if [ $# -eq 0 ]
then
     	usage
	exit 1
fi

# Parse options
FORCE=false

while getopts 'f' option
do
  	case $option in
                f) FORCE=true ;;
        esac
done

# skip to parse command line arguments
shift $((OPTIND-1))

SPLIT_BREAK_DIR=$(readlink -e $1) # absolute path to the split-break directory
PAIR_ID=$(basename $SPLIT_BREAK_DIR) # patient episode id

INPUT=$(readlink -e $SPLIT_BREAK_DIR/$PAIR_ID.mask.insertion.pairs.bed) # input file
# check that the input file is not empty. If true, stop

if [[ ! -s $INPUT ]]; then
    echo "The input file $INPUT is empty. Exiting run-annotate.splitters.sh" 
    exit 1
fi

# check if directory exist. If true and --force is selected, remove
if [ -d $PAIR_ID ]
then
	if $FORCE
		then echo "Removing existing $PAIR_ID directory"
                rm -rf $PAIR_ID
        else
            	echo "Directory $PAIR_ID exists. Use -f to remove directory. Exiting"
                exit 1
	fi
fi

# create file structure
# pair directory
mkdir $PAIR_ID
cd $PAIR_ID
mkdir reference
cd reference
cp $SPLIT_BREAK_DIR/reference/ref.* .
genbank2gff.pl < ref.gbk > ref.gff
cd ..
cut -f1 $SPLIT_BREAK_DIR/isolates.tab > isolates.txt
cp $INPUT .

# isolate directories
for ISO in $(cut -f4 $INPUT | sort -u)
do
	mkdir $ISO
	cd $ISO
	echo "Extracting primary splitters intervals for $ISO"
	grep $ISO $INPUT > "${ISO}.insertion.pairs.bed"
	bedtools merge -i $SPLIT_BREAK_DIR/$ISO/$ISO.splitters.bam | bedtools coverage -a stdin -b $SPLIT_BREAK_DIR/$ISO/$ISO.splitters.bam | cut -f1,2,3,4,5 > $ISO.splitters.bed
	echo "Primary splitters intervals extracted"
	cd ..
done


# run python script
cd ..
eval "$(conda shell.bash hook)"
conda activate base

python ~/perl5/bin/annotate-splitters_v3.py -o $PAIR_ID -i $INPUT -d $SPLIT_BREAK_DIR
