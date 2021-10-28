
USAGE="Error: $0 takes 3 arguments \\n
Usage: $0 <query> <database> <out>"

if [ $# -lt 3 ]
then
	echo -e $USAGE
	exit 1
fi

# Get command line argument
QUERY=$1
DATABASE=$2
OUT=$3
N_THREADS=8

# Functions
exe() { echo "\$ $@" ; "$@" ; }

# Database
echo "Creating blast database for $DATABASE"
makeblastdb -dbtype prot -in $DATABASE

# Run blastp with default output format (for inspection of the local alignments)
echo "Running:"
exe blastp -query $QUERY -db $DATABASE -evalue 1e-5 -num_threads $N_THREADS -out $OUT.txt

# Run blastp with a tabular output format (for R)
FMT="6 std qlen"
echo "Running:"
exe eval blastp -query $QUERY -db $DATABASE -evalue 1e-5 -num_threads $N_THREADS -out $OUT.tab -outfmt \"$FMT\"

# Delete database files
rm $DATABASE.*

# Final statements
echo "$0 done"
echo "You can find the output in $OUT.txt and $OUT.tab"
