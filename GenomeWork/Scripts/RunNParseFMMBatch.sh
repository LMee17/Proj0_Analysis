#!usr/bin/bash

#USAGE bash RunFMMBatch.sh <psglist> <Alignment folder> < species list> <directory of paml predicted tree files>

#PSGlist = list of genes under positive selection
#Alignment Folder = folder of gene alignments fed into paml
#species list = A text file with the list of species included in the analyses. Make sure naming matches that used in the alignment
#Make sure that all directories are written in full path format, ie Alignments as Alignments/

mkdir ResultFiles/

#collect tree
while read p; do
	echo $p
	#find a tree file
	cp $5$p*tree .
	#remove branch designations
	sed -i 's/$1//g' $p*tree
	#Now for the alignments
	#cp it in
	cp $2$p*paml .
	#Need to convert the format from paml to FMM readable
	#remove the top line
	sed -i '1d' $p*paml
	#add a carrot before all species
	while read s; do
		sed -i "s/$s/\>$s/g" $p*paml
	done < $3
	#run fmm
	hyphy fmm --alignment $p*paml --tree $p*tree --output ResultFiles/$p.FMM.json
	rm *tmp
	rm $p*tree
	rm $p*paml
done < $1

#Make a table with gene followed by whether or not it is considered better to use FMM

printf "%s\t%s\t%s\t%s\n" "Gene" "Double-hit_vs_single-hit" "Triple-hit_vs_double-hit" "Triple-hit_vs_single-hit" > FMMResults.tsv

#extract results

for j in ResultFiles/*json; do
        gene=$(echo $j | awk -F "/" '{print$NF}' | awk -F "." '{print$1}')
        echo $gene
        dvs=$(grep "p-value" $j | head -n 1 | awk -F ":" '{print$2}')
        tvd=$(grep "p-value" $j | head -n 2 | sed '1d' | awk -F ":" '{print$2}')
        tvs=$(grep "p-value" $j | head -n 3 | sed '1,2d' | awk -F ":" '{print$2}')
        printf "%s\t%s\t%s\t%s\n" $gene $dvs $tvd $tvs >> FMMResults.tsv
done

#extract evidence ratios per site
#sites above 5 are considered CMD


for j in ResultFiles/*json; do
        gene=$(echo $j | awk -F "/" '{print$NF}' | awk -F "." '{print$1}')
        grep -A 1 "Three-hit" $j | sed 's/]//g' | sed 's/\[//g' | sed '1d' | sed 's/, /\n/g' > $gene.3hit.tmp
        count=1
        printf "%s\t%s\n" "Codon" "EvidenceRatio" > $gene.3hit.Results.tsv
        while read e; do
                printf "%s\t%s\n" "$count" "$e" >> $gene.3hit.Results.tsv
                let count=count+1
        done < $gene.3hit.tmp
        grep -A 1 "Two-hit" $j | sed 's/]//g' | sed 's/\[//g' | sed '1d' | sed 's/, /\n/g' > $gene.2hit.tmp
        count=1
        printf "%s\t%s\n" "Codon" "EvidenceRatio" > $gene.2hit.Results.tsv
        while read e; do
                printf "%s\t%s\n" "$count" "$e" >> $gene.2hit.Results.tsv
                let count=count+1
        done < $gene.2hit.tmp
        rm *tmp
done

mkdir ParsedResults/
mv *tsv ParsedResults/

