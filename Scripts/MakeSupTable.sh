#!usr/bin/bash

#Parse the complete results file (CodeML_Results_All_Nov21.tsv) into three tables - Canon / Non-Canon / Background
#Row names = Genes, columns = tested lineages, values = adjpvalue of LRT 

file=$1

#set up our lists to iterate through
awk -F "\t" '{print$7}' $file | sed '1d' | sort -u > class.tmp
awk -F "\t" '{print$5}' $file | sed '1d' | sort -u > lineage.tmp

#per class, run through the file and pull out the results for that class of genes
while read c; do
	echo $c
	#Remove any spaces for file name
	cfile=$(echo $c | sed -e 's/\ /\_/g')
	echo $cfile
	#make a file ready to populate
	out="$cfile.BranchSite.CodeMLData.Nov21.tsv"
	printf "#Each branch site analysis results in a LRT score that is compared to a chi-square distribution with 1 degree of freedom to give a pvalue\n" >> $out
	printf "#This table shows the adjusted pvalues for each lineage per tested gene. NA means the test resulted in convergence failure\n\n" >> $out
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Gene" "All_Post_Elaboration_Sociality" "All_Post_Origin_Sociality" "All_Solitary" "Apis" "Social_Corbiculates" "Habropoda" "Lasioglossum" "Megachilidae" "Melipona" "Dufourea" "Ceratina" >> $out
	grep "\t$c\t" $file > $cfile.tmp
	#now we go through by gene
	awk -F "\t" '{print$1}' $cfile.tmp | sed '1d' | sort -u > $cfile.gene.tmp
	tot=$(wc -l $cfile.gene.tmp | awk '{print$1}')
	count=1
	while read g; do
		echo "$g, $count / $tot"
		#separate out each individual gene results from the class file
		grep $g	$cfile.tmp > $cfile.$g.tmp
		#and then go through by lineage
		while read l; do
			#check the sociality is there
			check=$(grep -c "$l\t" $cfile.$g.tmp)
			if [ $check -eq 1 ]; then
				p=$(grep "$l\t" $cfile.$g.tmp | awk -F "\t" '{print$15}')
				printf "%s\n" "$p" >> $cfile.$g.pvalue.tmp
			else
				printf "%s\n" "NA" >> $cfile.$g.pvalue.tmp
			fi 
		done < lineage.tmp
		printf "%s\t" "$g" >> $out 
		while read r; do
			printf "%s\t" "$r" >> $out
		done < $cfile.$g.pvalue.tmp 
		printf "\n" >> $out
		rm $cfile.$g.pvalue.tmp
		rm $cfile.$g.tmp
		let count=count+1
	done < $cfile.gene.tmp
done < class.tmp

rm *tmp
