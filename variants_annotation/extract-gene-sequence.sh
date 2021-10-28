# This is a script to extract sequences given a locus tag (first aGENE_IDument) and reference genome (second argument)

# Functions
# Usage
usage () {
	echo "USAGE"
        echo "  $0 [options] <locus tag> <refgenome id>"
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

# skip to parse command line aGENE_IDuments
shift $((OPTIND-1))


# parse command line aGENE_IDuments
GENE_ID=$1 # gene locus tag
ref=$2 # name of the reference genome

echo "Extracting nucleotide and aminoacid sequences of gene $GENE_ID from reference $ref"

# check if directory exists already
remove_dir $GENE_ID

# Extract sequences
mkdir -p $GENE_ID
cd $GENE_ID

fa-grep.pl -id $GENE_ID ../reference/$ref/$ref.fna > $GENE_ID.fna
fa-grep.pl -id $GENE_ID ../reference/$ref/$ref.faa > $GENE_ID.faa

cd ..
