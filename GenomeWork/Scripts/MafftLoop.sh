#!usr/bin/bash

# Loop through a number of multiple sequence fasta files in a directory and create alignments 

# USAGE bash MafftLoop.sh <directory containing fastas> <fasta extension> <nuc/amino>

for f in $1/*$2; do
	ID=$(basename "$f" .$2)
	echo "$ID"
	mafft --globalpair --maxiterate 1000 --$3 --treeout --clustalout $f > $ID.aln
done

echo "Alignment runs complete".
