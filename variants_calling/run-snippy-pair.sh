# Assess arguments
if [ $# -lt  3 ]
then
	echo "USAGE:"
	echo " $0 [options] <PAIR ID> <iso1> <iso2>"
	echo "OPTIONS:"
	echo " -f: force"
	echo " -T(threads): (default =8)"
	echo " -r: path to reference draft assemblies"
	echo " -i: path to the isolates.tab file with the reads paths"
	exit 1
fi

# Parse options
FORCE=false
THREADS=8

while getopts 'fT:r:i:' option
do
	case $option in
		f) FORCE=true ;;
		T) THREADS=$OPTARG ;;
		r) REFPATH=$OPTARG ;;
		i) ISOFILE=$OPTARG ;;
	esac
done

# skip to parse command line arguments
shift $((OPTIND-1))


# Get variables
PAIR_ID=$1
iso1=$2
iso2=$3
ISOFILE=$(readlink -e $ISOFILE)

# Starting statements
echo "Running $0 on $PAIR_ID"
echo "Reference (iso1) is $iso1"
echo "iso2 is $iso2"

# Create pair directory
if [ -d $PAIR_ID ]
then
	if $FORCE
	then
		echo "$PAIR_ID exists. $0 will remove $PAIR_ID"
		rm -r $PAIR_ID
	else
		echo "$PAIR_ID exists. Use --f to remove directory. Exiting"
		exit 1
	fi
fi

mkdir $PAIR_ID
cd $PAIR_ID

 # Create isolates.tab file
grep $iso1 $ISOFILE >> isolates.tab
grep $iso2 $ISOFILE >> isolates.tab

# Copy reference
mkdir reference 
cp $REFPATH/$iso1/$iso1.gbk reference

# Run snippy
cat isolates.tab | while read f1 f2 f3;
	do echo "Running snippy on $f1"
	snippy --cpus $THREADS --ref reference/*.gbk --outdir $f1 --R1 $f2 --R2 $f3
done

cd ..

# Final statements
echo "$0 on $PAIR_ID completed"
echo "You can run-concatenate-snippy.sh to remove variants from the reference reads"
echo "Or, even better, you can run mask-variants.sh to remove variants from the reference reads and variants in positions where reference reads have low coverage"
echo "This will substantially improve the accuracy of your mutation list!"
