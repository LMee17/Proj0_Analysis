#!usr/bin/bash

# USAGE bash Pal2NalLoop.sh <directory containing aln> <directory containing fasta> <output folder>

# SCRIPT TO ITERATE THROUGH A DIRECTORY OF PROTEIN ALIGNMENTS AND A CORRESPONDING DIRECTORY OF MULTIPLE TRANSCRIPT FASTA FILES
# IN ORDER TO CREATE CODON ALIGNMENTS READY TO BE USED IN PAML. REQUIRES ALL PROTEIN ALIGNMENTS TO END IN *ALN AND ALL FASTA FILES
# TO END IN *FASTA. ASSUMES THAT EACH MATCHING PEPTIDE AND FASTA FILE HAS THE SAME NAME BEFORE THE FILE EXTENSION. THIS SCRIPT 
# WILL ALSO SEPARATE ANY FAILED ALIGNMENTS INTO ANOTHER FOLDER WITHIN THE OUTPUT CALLED FAILED - A .TXT FILE OF ALL FAILED
# ALIGNMENTS WILL ALSO BE CREATED FOR EASE OF USE IN DOWNSTREAM PROCESSES (IE BLAST RE-RUNS OR DATASET CULLS)

mkdir $3/

for a in $1/*.aln; do
	ID=$(basename "$a" .aln)
	echo "$ID"
	perl ~/bin/pal2nal.v14/pal2nal.pl $1/$ID.aln $2/$ID.fasta -output paml -nogap > $ID.paml
	COUNT=$(wc -l $ID.paml | awk '{print $1}' )
	echo $COUNT
	if [ $COUNT -eq 0 ]; then
		printf "%s\n" "$ID" >> pal2nalfail.txt
	else
		mv $ID.paml $3/
	fi
done

if [ -e pal2nalfail.txt ]; then
	mkdir $3/Failed/
	mv *paml $3/Failed/
	mv pal2nalfail.txt $3/Failed/
else
	echo "No Codon Alignment Failures"
fi

echo "Pal2Nal Codon Alignment Generation Complete"
