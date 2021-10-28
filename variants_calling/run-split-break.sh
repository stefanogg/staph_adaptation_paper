# Shell wrapper to run the split-break script on genetic pairs splitters bam

# It takes one argument: the directory with the splitters 

# Usage
usage () {
	echo "USAGE"
	echo "  $0 [options] <pair splitters dir> <internal reference>"
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


# Get arguments

SPLITTERS_DIR=$(readlink -e $1) # absolute path to the patient episode splitters directory
PAIR_ID=$(basename $SPLITTERS_DIR) # patient episode id

INPUT=$(readlink -e $SPLITTERS_DIR/splitters.tab)
REF=$(readlink -e $SPLITTERS_DIR/reference/*.gbk)
INTERNAL_REF=$2
#c=$(grep $p internal_ref.tab | cut -f2)

# Check if directory exists and delete if true
if [ -d $PAIR_ID ]
then
	if $FORCE
	then 
		echo "Removing existing $PAIR_ID directory"
		rm -rf $PAIR_ID
	else
		echo "Directory $PAIR_ID exists. Use -f to remove directory. Exiting"
		exit 1
	fi
fi

# Run script
python ~/perl5/bin/split-break_v2.py -o $PAIR_ID -i $INPUT -r $REF -c $INTERNAL_REF

# Generate input file for annotate-splitters
cd $PAIR_ID

for i in $(cut -f1 isolates.tab)
	do cut -f1,2,10 $i/$i.mask.insertion.pairs.bed | sed "s/$/\t$i/" > $PAIR_ID.mask.insertion.pairs.bed
done

# Mask genetically unrelated strains
#for i in $(grep $PAIR_ID ../close_related.tab | cut -f2)
#	do cut -f1,2,10 $i/$i.mask.insertion.pairs.bed | sed "s/$/\t$i/" > $p.close.mask.insertion.pairs.bed
#done 

cd ..
