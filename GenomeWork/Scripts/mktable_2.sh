#!usr/bin/bash

#USAGE bash <script> <directory of codon alignments *paml> 


#must be ran in same directory as BranchSiteCodeml.sh result file [name of run].report with no other files
#with the same suffix

wd=$(pwd | awk -F "/" '{print$(NF-1)}')
echo $wd
run=$(pwd | awk -F "/" '{print$NF}')
echo $run

printf "%s\t%s\t%s\n" "Gene" "InitialKappa" "EstTree" >> $run.M0.readout


for g in $1/*paml; do
	echo $g
	gene=$(echo $g | awk -F "/" '{print$NF}' | awk -F "." '{print$1}')
	echo $gene
	grep -A 3 "$gene" *.report > $gene.tmp
	kappa=$(grep "Kappa" $gene.tmp |  awk -F " " '{print $3}')
	echo $kappa
	tree=$(grep "):" $wd/$gene.$wd/$gene.$wd.M0.out | tail -n 1)
	echo $tree
	printf "%s\t%s\t%s\n" "$gene" "$kappa" "$tree" >> $run.M0.readout
	rm *tmp 
done
