#usr/bin/bash

#USAGE bash <script> <goi> <proteome> <table of gene/transcript/protein IDs> <isoform counts> <peptide name/gene ID file>  <output name>

# NEED TO SPLIT GENES INTO THOSE THAT ONLY HAVE ONE ISOFORM AND THOSE THAT HAVE MORE THAN ONE

sort -u $1 > inputgoi
printf "%s\n" "Beginning Resolution" >> $6.IsoResolve.report


while read g; do
	printf "%s\n" "$g" >> $6.IsoResolve.report
	count=$(grep -P "$g\t" $4 | awk '{print $2}')		#take the isoform count from a file provided specifying the number of protein isoforms per gene (produced by GFFread.sh)
	if [ $count == 0 ]; then 							# if there are no protein isoforms, assume that the gene does not produce any and skip.
		echo "$g has no protein product. Skipping."
		printf "%s\n" "No protein product found. Skipping." >> $6.IsoResolve.report
	else
		# IF THERE IS MORE THAN ONE ISOFORM FOR THIS GENE, THIS NEEDS TO BE RESOLVED.
		if [ $count -gt 1 ] ; then 					#THIS PART OF THE SCRIPT IS ONLY ACTIVATED IF THERE ARE MORE THAN ISOFORM
			echo "$g has more than one isoform. Resolving."
			printf "%s\n" "More than one isoform detected. Resolving" >> $6.IsoResolve.report
			grep "\<$g\>" $3 | awk '{print $3}' | sort -u | sort > $g.tmp 			# Extract the gene from the ID table provided | Extract just the protein product IDs | Check there are no duplicates | Sort alphabetically | Store in a temporary list file
			nccheck=$(grep -c "P" $g.tmp) 				# Count the number of protein products in the list. This is a failsafe check to ensure no noncoding RNA is extracted by accident.
			if [ $nccheck == 0 ]; then
				echo "$g does not code for a protein product. Skipping"
				printf "%s\n" "$g" >> $6.NonCodingGenes_Skipped.txt
				printf "%s\n" "No protein product for this gene could be detected with the given input files. Skipping." >> $6.IsoResolve.report
			else
				# IF THERE'S A CURATED ISOFORM (NP), THIS SHOULD BE PREFERRED OVER ALL OTHERS
				npcheck=$(grep -c "NP_" $g.tmp) 		#Check for presence of NP_ proteins
				# IF THERE'S NO CURATED ISOFORMS
				if [ $npcheck == 0 ]; then 				# If there are no curated isoforms present, then move on to using IsoSel to resolve amongst the multiple non-curated isoforms.
					echo "$g has no curated isoform present. Using IsoSel to resolve isoforms."
					printf "%s\n" "No curated isoforms present. Using IsoSel to resolve isoforms." >> $6.IsoResolve.report
					python ~/Scripts/get_seq2.py $2 $g.tmp $g.fas 			# Make a fasta file of the gene's multiple non-curated isoforms
					grep "\<$g\>" $5 > $g.iso.tmp 				# Make a locus/isoform tag file just for that gene using the master file provided above (produced by GFFread.sh)
					/usr/local/bin/run_isosel $g.fas -n $g.iso -auto -f $g.iso.tmp 			#Run isosel. The presence of the locus/isoform file makes IsoSel produce a filtered.fasta file with only the resolved isoform included
					failcheck=$(grep -c "End of the IsoSel run" $g.iso.log) 		#If IsoSel is successful, then this line is present in the log file it produces.
					# IF ISOSEL FAILS
					if [ $failcheck == 0 ]; then 			#If this line is not present, assume it has failed and move on to selecting by length instead.
						echo "IsoSel failed to resolve $g isoforms. Selecting longest isoform instead."
						while read p; do 		#Run through the temporary per gene protein list produced above.
							echo "$p" > $p.id.tmp 		# For each gene, make a temporary file with just the ID included
							python ~/Scripts/get_seq2.py $2 $p.id.tmp $p.fas.tmp 	#Use this single ID to produce fasta files with just the one isoform included
							sed -i '1d' $p.fas.tmp 				#Remove the first line from the fasta file, leading just protein sequence
							charcount=$(wc $p.fas.tmp | awk '{print $3}') 			#Count the characters within this temporary file, therefore counting the number of amino acids in the sequence.
							printf "%s\t%d\n" "$p" "$charcount" >> $g.comp.tmp 		#This character count is stored and then printed into comparative file with peptide \t sequence length
						done < $g.tmp
						sort -nk 2 $g.comp.tmp | tail -n 1 | awk '{print $1}' > $g.longiso.tmp 		# This comparative file is then sorted so that the longest isoform length is bottom. -n = sort by number, -k = sort by a particular column (2 in this case). Take the bottom most line (the longest) and extract that protein ID
						python ~/Scripts/get_seq2.py $2 $g.longiso.tmp $g.iso_filtered.fasta 	#Make a fasta file with just this longest isoform
						rm $g.comp.tmp 			#Remove the comparative file and anything else produced by IsoSel not needed for the next stage
						rm *fas
						rm *log
						rm *aln
						printf "%s\n" "IsoSel failed to resolve isoforms. Selecting longest isoform." >> $6.IsoResolve.report
					# IF ISOSEL WORKS
					else
						echo "IsoSel successfully resolved $g's isoforms."
						printf "%s\n" "IsoSel resolved isoforms." >> $6.IsoResolve.report
						rm *log 			#Clean up
						rm *Scores
						rm *aln
						rm *scores
					fi
				else 
					# IF THERE'S JUST ONE CURATED ISOFORM
					if [ $npcheck == 1 ]; then 
						echo "Selecting $g's curated isoform."
                        			grep "NP_" $g.tmp > $g.iso.tmp
                        			python ~/Scripts/get_seq2.py $2 $g.iso.tmp $g.iso_filtered.fasta # Just make a fasta of that one NP - this is considered the iso_filtered.fasta
						printf "%s\n" "Curated isoform present. Selecting." >> $6.IsoResolve.report
					# IF THERE'S MORE THAN ONE CURATED ISOFORM
					else
						echo "Resolving between $g's curated isoforms"
						printf "%s\n" "More than one curated isoform detected. Using IsoSel to resolve" >> $6.IsoResolve.report
						grep "NP_" $g.tmp > $g.NP.tmp 		#Pull out all the curated isoforms and put them in a temporary list file
						grep "\<$g\>" $5 | grep "NP_" > $g.iso.tmp 		#Pull out the isoform / gene locus information for these isoforms
						python ~/Scripts/get_seq2.py $2 $g.NP.tmp $g.fas 		#Make a fasta file of these curated isoforms
						/usr/local/bin/run_isosel $g.fas -n $g.iso -auto -f $g.iso.tmp 		#Use the isoform/locus and fasta file information to get isosel to resolve isoforms
						failcheck=$(grep -c "End of the IsoSel run" $g.iso.log) 		#Check to see if IsoSel ran to completion
						rm *fas
						# IF SISOSEL FAILS
		        	       		if [ $failcheck == 0 ]; then 		#If this line is not present, assume it has failed and move on to selecting by length instead.
                					echo "IsoSel failed to resolve between curated $g isoforms. Selecting longest isoform instead."
                        				while read p; do 		#Run through the temporary per gene protein list produced above.
                               					echo "$p" > $p.id.tmp 		# For each gene, make a temporary file with just the ID included
                                				python ~/Scripts/get_seq2.py $2 $p.id.tmp $p.fas.tmp 		#Use this single ID to produce fasta files with just the one isoform included
                                				sed -i '1d' $p.fas.tmp 			#Remove the first line from the fasta file, leading just protein sequence
                                				charcount=$(wc $p.fas.tmp | awk '{print $3}') 		#Count the characters within this temporary file, therefore counting the number of amino acids in the sequence.
                                				printf "%s\t%d\n" "$p" "$charcount" >> $g.comp.tmp 			#This character count is stored and then printed into comparative file with peptide \t sequence length
                            				done < $g.NP.tmp
                            				sort -nk 2 $g.comp.tmp | tail -n 1 | awk '{print $1}' > $g.longiso.tmp 		# This comparative file is then sorted so that the longest isoform length is bottom. -n = sort by number, -k = sort by a particular column (2 in this case). Take the bottom most line (the longest) and extract that protein ID
                            				python ~/Scripts/get_seq2.py $2 $g.longiso.tmp $g.iso_filtered.fasta 		#Make a fasta file with just this longest isoform
							rm $g.comp.tmp 			#Remove the comparative file and anything else produced by IsoSel not needed for the next stage
							rm *log 			#Cleanup
							rm *aln
							printf "%s\n" "IsoSel failed to resolve between curated isoforms. Choosing longest curated isoform." >> $6.IsoResolve.report
						else 
							echo "IsoSel successfully resolved between $g curated isoforms"
							rm *log 		#Cleanup
							rm *Scores
							rm *aln
							rm *scores
							printf "%s\n" "IsoSel resolved between curated isoforms." >> $6.IsoResolve.report
						fi
					fi
				fi
           		 fi
		else
			# IF THERE'S ONLY ONE NON CURATED ISOFORM
			echo "$g does not have multiple isoforms"
			grep "\<$g\>" $3 | awk '{print $3}' | sort -u | sort > $g.tmp  #Make a temporary file with the id of this isoform
			python ~/Scripts/get_seq2.py $2 $g.tmp $g.iso_filtered.fasta 		#Use this id to make a fasta file - this will be added as is to the final filtered.faa
			printf "%s\n" "There are no multiple isoforms to resolve. Proceeding with single protein product." >> $6.IsoResolve.report
		fi
		rm *tmp
	fi
done < inputgoi 		#All of this takes place whilst the script iterates through the genes in the gene list provided.

rm inputgoi

printf "%s\t%s\n" "Gene" "Resolved_Isoform" >> $6.Filteredfaa.table 		#Make a table for easy access of gene / resolved isoform
for i in *iso_filtered.fasta; do
	iso=$(grep ">" $i | awk -F ">" '{print $2}' | awk '{print $1}') 		#Go through each filtered isoform and record the protein ID 
	gene=$(basename $i ".iso_filtered.fasta") 								#Record the name of the gene from the file name.
	printf "%s\t%s\n" "$gene" "$iso" >> $6.Filteredfaa.table 			#Put into a table
done

cat *_filtered.fasta > Isoform_filtered.tmp
#mkdir FilteredFastas/
#mv *_filtered.fasta FilteredFastas/

grep ">" Isoform_filtered.tmp | awk -F ">" '{print $2}' | awk '{print $1}' | sort -u > protlist.tmp
python ~/Scripts/get_seq2.py $2 protlist.tmp "$6_Isoform_filtered.faa"

rm *tmp
rm *_filtered.fasta

#mkdir FilteredFastas/
#mv *_filtered.fasta FilteredFastas/

echo "Isoform Resolution Complete"
