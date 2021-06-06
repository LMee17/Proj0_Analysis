#!usr/bin/bash

#Parse the complete results file (CodeML_BranchSiteAnalysis_CompleteResults.tsv) into three tables - Canon / Non-Canon / Background
#Row names = Genes, columns = tested lineages, values = adjpvalue of LRT 

#combine the two results files
cp CodeML_BranchSiteAnalysis_Sol_CompleteResults.tsv solresults.tmp
sed -i '' '1d' solresults.tmp
cat CodeML_BranchSiteAnalysis_CompleteResults.tsv solresults.tmp > CodeML_BranchSiteAnalysis_EntireData.tsv

rm solresults.tmp

#that's a typeful so
file=CodeML_BranchSiteAnalysis_EntireData.tsv

#make sure NonCanon is changed to Non-Canon (as is used in the manuscript)
sed -i '' 's/NonCanon/Non-Canon/g' $file

#set up our lists to iterate through
awk -F "\t" '{print$6}' $file | sed '1d' | sort -u > class.tmp
awk -F "\t" '{print$5}' $file | sed '1d' | sort -u > lineage.tmp

#put out the silly episodic runs that aren't being used
sed -i '' '/AllOriginEpi/d' lineage.tmp
sed -i '' '/CorbSocEpi/d' lineage.tmp


#per class, run through the file and pull out the results for that class of genes
while read c; do
	echo $c
	#make a file ready to populate
	printf "#Each branch site analysis results in a LRT score that is compared to a chi-square distrubtion with 1 degree of freedom to give a pvalue\n" >> $c.BranchSite.CodeMLData.tsv
	printf "#This table shows the adjusted pvalues for each lineage per tested gene. NA means the test resulted in convergence failure\n\n" >> $c.BranchSite.CodeMLData.tsv
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Gene" "All_Advanced_Eusocial" "All_Social" "All_Solitary" "Apis" "Social_Corbiculates" "Habropoda" "Lasioglossum" "Megachilidae" "Melipona" "Dufourea" "Ceratina" >> $c.BranchSite.CodeMLData.tsv
	out="$c.BranchSite.CodeMLData.tsv"
	grep "\t$c\t" $file > $c.tmp
	#now we go through by gene
	awk -F "\t" '{print$1}' $c.tmp | sed '1d' | sort -u > $c.gene.tmp
	tot=$(wc -l $c.gene.tmp | awk '{print$1}')
	count=1
	while read g; do
		echo "$g, $count / $tot"
		#separate out each individual gene results from the class file
		grep $g	$c.tmp > $c.$g.tmp
		#and then go through by lineage
		while read l; do
			#check the sociality is there
			check=$(grep -c "$l\t" $c.$g.tmp)
			if [ $check -eq 1 ]; then
				p=$(grep "$l\t" $c.$g.tmp | awk -F "\t" '{print$9}')
				printf "%s\n" "$p" >> $c.$g.pvalue.tmp
			else
				printf "%s\n" "NA" >> $c.$g.pvalue.tmp
			fi 
		done < lineage.tmp
		printf "%s\t" "$g" >> $out 
		while read r; do
			printf "%s\t" "$r" >> $out
		done < $c.$g.pvalue.tmp 
		printf "\n" >> $out
		rm $c.$g.pvalue.tmp
		rm $c.$g.tmp
		let count=count+1
	done < $c.gene.tmp
done < class.tmp

rm *tmp
