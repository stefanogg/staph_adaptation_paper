# Shell pipeline to map reads on the final ref genome using bwa mem

# Usage

USAGE="$0 [option] <pair id> <refgenome> <iso1> <iso2>\\n
-f (force): remove existing directory\\n
-t (threads): number of threads (for bwa) (default: 8)"

if [ $# -eq 0 ]
then
	echo -e $USAGE
	exit 1
fi

# initial statements
echo "This is map_on_ref_genomes.sh"
echo "This script uses bwa to map reads to the final reference genome"

# parse options
FORCE=false
THREADS=8

while getopts 'ft:' option
do
	case $option in 
		f) FORCE=true ;;
		t) THREADS=$OPTARG ;;
	esac
done

echo "Number of threads (for bwa and samtools): $THREADS"

# skip to parse command line arguments
shift $((OPTIND-1))

# Parse command line arguments
PAIR_ID=$1
REFGENOME=/home/giulieris/saureus-general/ref-genomes/download_mar21/${2}
ISO1=$3
ISO2=$4
ISOFILE=$(readlink -e isolates.tab)
LOG="map_on_ref_genomes.log"

# Create file structure

# check if group directory exists and remove if true (later to be replaced by --force option)
if [ -d $PAIR_ID ]
then
       	if $FORCE
       	then
               	echo "Directory $PAIR_ID exists. Removing"
               	rm -r $PAIR_ID
       	else
               	echo "Directory $PAIR_ID exists. Use -f to remove directory. Exiting"
               	exit 1
       	fi
fi

# create directory
mkdir -p $PAIR_ID
cd $PAIR_ID

echo "Running map_on_ref-genomes.sh on pair $PAIR_ID" 2>&1 | tee -a $LOG


# copy relevant reference files
REFNAME=$(basename -s .fna.gz $REFGENOME)
REFDIR=$(dirname $REFGENOME)
echo "Copying reference $REFNAME genbank file" 2>&1 | tee -a $LOG

mkdir reference
cd reference
# copy refgenome in genbak format
cp ${REFDIR}/$REFNAME.gbff.gz .
gunzip *
mv $REFNAME.gbff $REFNAME.gbk
any2fasta -u $REFNAME.gbk > $REFNAME.fna

cd ..


# Create isolates.tab file
grep $ISO1 $ISOFILE >> isolates.tab
grep $ISO2 $ISOFILE >> isolates.tab

# run bwa

cat isolates.tab | while read f1 f2;
	do echo "Mapping reads of strain $f1"
	mkdir $f1
	cd $f1
	mkdir reference
	cd reference
	cp ../../reference/*fna ref.fa
	bwa index ref.fa
	cd ..
	bwa mem -t $THREADS reference/ref.fa $f2 | samtools sort -@ $THREADS > $f1.bam
	samtools index $f1.bam
	cd ..
done 2>&1 | tee -a $LOG

# Final statements
echo "[$(date)] $0 done" | tee -a $LOG

cd ..

