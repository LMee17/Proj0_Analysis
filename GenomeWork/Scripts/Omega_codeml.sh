#!usr/bin/bash

#USAGE: bash Omega_codeml.sh <directory containing alignments> <ctl file> <table with M0 estimates>


printf "%s\t%s\n" "Gene" "dN/dS" >> dNdSratio.tsv
mkdir RunFiles/

check=$(echo $1 | grep -c "/")
echo $1
if [ $check -eq 1 ]; then
        WD=$(echo $1 | sed 's/.$//')
        echo $WD
else
        WD=$1
        echo $WD
fi

WD=$(echo "$WD" | sed 's/\//\\\//g')
echo $WD

for p in $1/*paml; do
	#set up ctl file per gene and run codeml
	gene=$(echo $p | awk -F  "/" '{print$NF}' | awk -F "." '{print$1}')
	file=$(echo $p | awk -F "/" '{print$NF}')
	echo $gene
	echo $file
	sed -i "s/seqfile = .*/seqfile = $WD\/$file/" $2
	sed -i "s/outfile = .*/outfile = $gene.omega.out/" $2
	sed -i "s/treefile = .*/treefile = $gene.tree/" $2
	grep $gene $3 | awk -F "\t" '{print $3}'> $gene.tree
	kappa=$(grep $gene $3 | awk -F "\t" '{print $2}')
	echo $gene
	sed -i "s/ kappa = .*/ kappa = $kappa/" $2
	codeml $2
	#extract what we need
	omega=$(grep "(dN/dS)" $gene.omega.out | awk -F "=" '{print $2}')
	#record
	printf "%s\t%s\n" "$gene" "$omega" >> dNdSratio.tsv
	#Cleanup
	mv 2NG.dN $gene.omega.2NG.dN
	mv 2NG.dS $gene.omega.2NG.dS
	mv 2NG.t $gene.omega.2NG.t
	mv lnf $gene.omega.lnf
	mv rst $gene.omega.rst
	mv rst1 $gene.omega.rst1
	mv rub $gene.omega.rub
	mkdir RunFiles/$gene/
	mv $gene* RunFiles/$gene/
done
