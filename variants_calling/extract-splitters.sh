# Usage

USAGE="$0 <grouped bwa directory>"

if [ $# -eq 0 ]
then
	echo $USAGE
	exit 1
fi


BWA_DIR=$(readlink -e $1) # bwa directory with bam files
PAIR_ID=$(basename $BWA_DIR) # patient episode unique id
LOG="$PAIR_ID.extract-splitters.log"

[[ -d $PAIR_ID ]] && rm -r $PAIR_ID # remove directory if existing

mkdir $PAIR_ID
cd $PAIR_ID

echo "Running extract-splitters on $PAIR_ID" 2>&1 | tee -a $LOG

cp -r $BWA_DIR/reference .

for i in $(cut -f1 $BWA_DIR/isolates.tab)
	do echo "Extracting split reads of strain $i"
	mkdir $i
	cd $i
	samtools view -h $BWA_DIR/$i/$i.bam | python /home/giulieris/perl5/bin/extractSplitReads_BwaMem.py -i stdin | samtools view -Sb | samtools sort > $i.splitters.bam
	samtools index $i.splitters.bam
	r=$(readlink -e $i.splitters.bam)
	echo -e "$i\t$r" >> ../splitters.tab
	cd ..
done 2>&1 | tee -a $LOG

cd .

