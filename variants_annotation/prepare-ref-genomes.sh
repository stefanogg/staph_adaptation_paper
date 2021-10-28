# Script to extract faa and fna gene files from a list of reference genomes

# Functions
# Usage
usage () {
	echo "USAGE"
        echo "  $0 [options] <refgenome id> </path/to/refgenome>"
        echo "OPTIONS"
        echo "  -f (force): remove existing directory"
}

# Remove existing directory
remove_dir () {
	OUTDIR=$1
	if [ -d $OUTDIR ]
	then 
		if $FORCE
		then
			echo "Removing existing $OUTDIR directory"
			rm -rf $OUTDIR
		else
			echo "Directory $OUTDIR exists. Use -f to remove directory. Exiting $0"
			exit 1
		fi
	fi
}

# Create gene sequence files
prepare_refgenome () {
	cp $2 $1.gbk	
	genbank2fasta.pl $1.gbk > $1.fna
	genbank2fasta.pl --prot $1.gbk > $1.faa
}

# Run script

# Print usage
if [ $# -lt 2 ]
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

# Parse positional arguments
REFGENOME_ID=$1
REFGENOME_PATH=$(readlink -e $2)

# Create refgenome directory
remove_dir $REFGENOME_ID
mkdir -p $REFGENOME_ID

# Extract gene sequences in faa and fna format
cd $REFGENOME_ID
prepare_refgenome $REFGENOME_ID $REFGENOME_PATH

# Final statements
cd ..
echo "[$(date)] $0 done"
