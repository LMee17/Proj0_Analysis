#!usr/bin/bash

# USAGE bash <script> <list of transcripts> <anchor species> <no_threads> <working directory from ~/ *path to* > <no of species>
# REQUIREMENTS
# PAL2NAL AND SCRIPT Pal2NalLoop.sh (IN DIRECTORY ~/Scripts/)
# MAFFT AND SCRIPT MafftLoop.sh (IN DIRECTORY ~/Scripts/)
# getseq2.py SCRIPT (IN DIRECTORY ~/Scripts/)
# BLAST+
# MAKE A WORKING DIRECTORY. PROVIDE THE PATH ABOVE FROM ~ (IE ~/WORK/TEST/ = WORK/TEST)
# A DIRECTORY FULL OF PROTEOMES WITHIN THE WORKING DIRECTORY. 
# ALL PROTEOMES ARE ASSUMED TO HAVE SUFFIX *FAA. 
# A DIRECTORY FULL OF TRANSCRIPTOMES WITHIN THE WORKING DIRECTORY. 
# ALL TRANSCRIPTOMES ARE ASSUMED TO HAVE SUFFIX *FASTA
# ALSO REQUIRES A MASTER CDS // PROTEOME WITH SPECIES AS ID FOLLOWED BY ID IN RESPECTIVE # DIRECTORIES
# CDS NAME = MASTERCDS.FASTA, PEP NAME = MASTERPEP.FA
# THIS SCRIPT ALSO ASSUMES THAT THERE IS A MASTER PEP FILE WITH IDS STILL. IF THIS IS
# NOT THE CASE UNHASH THE LINE BELOW WHERE ONE IS MADE 
# ALL TRANSCRIPTOMES AND PROTEOMES MUST BE NAMED WITH A FOUR LETTER CODE AS PREFIX
# IE APIS MELLIFERA = AMEL_ETC
# THIS SCRIPT WILL NOT WORK IF RNA AND PROTEIN IDS ARE THE SAME. ONE WAY TO ENSURE THIS ISN'T THE CASE
# IS TO ENSURE THAT ONE (IE PROTEIN) IS CAPITALISED AND THE OTHER IN LOWER CASE BEFORE RUNNING THE SCRIPT
# THE LIST OF TRANSCRIPTS MUST BE FROM THE ANCHOR SPECIES. 
# ANCHOR SPECIES IDS ARE ASSUMED TO BE IN NCBI FORMAT, IE NM_0000154.2

file=$1


# MAKE FASTA FILE OF THE WANTED TRANSCRIPTS
echo "Producing fasta file from wanted transcripts"
cp ~/$4/Transcriptomes/$2*.fasta .                                          #Get transcriptomes into WD
python ~/Scripts/get_seq2.py $2*.fasta $1 wanted.fa                         #Make new fasta file of all wanted transcripts from anchor species transcriptome

IN=$(wc -l $1 | awk '{print $1}')                                           #Make a note of how many transcripts were first input
WANTNO=$(grep -c ">" wanted.fa)                                             #Make a note of how many were recovered from the anchor transcriptome

if [ $IN -eq $WANTNO ]; then
	echo "$WANTNO Transcripts successfully recovered"
else
	echo "$WANTNO Transcripts recovered"
	DIFF=$(($IN - $WANTNO))
	echo "$DIFF Transcripts were not recovered"
fi

# BLASTX TO GET PROTEIN PRODUCTS OF EACH TRANSCRIPT

echo "Blasting transcripts against anchor species proteome"
cp ~/$4/Proteomes/$2*.faa .                                                 #Get proteomes into WD
makeblastdb -in $2*.faa -out $2.prot.db -dbtype prot -parse_seqids          #Make BLAST databases for anchor species
blastx -query wanted.fa -db $2.prot.db -num_threads $3 -evalue 1e-10 -max_target_seqs 1 -outfmt 6 -out translate.out -best_hit_overhang 0.1 -best_hit_score_edge 0.1
awk '{print $2}' translate.out > translated.txt                             #Create a file of the translated transcripts' protein product
TRANSNO=$(wc -l translated.txt | awk '{print $1}')                          #Note how many were successfully translated

if [ $IN -eq $TRANSNO ]; then
	echo "$TRANSNO Transcripts successfully translated"
else
	echo "$TRANSNO Transcripts translated"
	DIFF=$(($IN - $TRANSNO))
	echo "$DIFF Transcipts failed to be translated"
fi

while read t; do                                                            #Go through the list of target transcripts line by line, pulling out the transcript and its translated product
        grep "$t" translate.out | awk '{print $2}' > $t.pep                 #from the blastx results and creating a new temp file named after the transcript, containing the peptide as the first line.
done < $file

# MAKE FASTA FILE OF PROTEINS

echo "Generating protein fasta"

python ~/Scripts/get_seq2.py $2*.faa translated.txt protein.fa              #Create a protein fasta of the list of translated products


# BLAST AGAINST ALL OTHER SPECIES

echo "Blasting protein fasta against other species"

rm $2*.faa                                                                 
rm $2*.fasta

cp ~/$4/Proteomes/*.faa .                                                   #Bring in all species' proteomes to begin RBH of anchor species peptides vs all other species.
cp ~/$4/Proteomes/*.fa .

cat *.faa > MasterPepID.fa                                                  #Produce a MasterPeptide file that contains all species and therefore all IDs

for f in *.faa; do                                                          #Make blast databases of all species
        echo $f
        makeblastdb -in $f -out $f.db -dbtype prot -parse_seqids
done

for b in *.faa; do                                                          #Blast anchor species' wanted peptides against other species (AvB run)
        echo $b
        blastp -query protein.fa -db $b.db -out $b.blastout -outfmt 6 -max_target_seqs 1  -num_threads $3 -evalue 1e-6 -best_hit_overhang 0.1 -best_hit_score_edge 0.1
done

cat *.blastout > protein.AvB.out                                            #Combine results
rm *.blastout                                                               #Clean up unnecessary files

echo "Blasting potential orthologs against anchor species"

awk '{print $2}' protein.AvB.out > Btarget.txt                              #Intermediate file consisting of all best hits from non-anchor species

for f in *.faa; do                                                          #Run through all species' proteomes and pull out best hit peptides to put back into a fasta file to blast back against anchor 
        echo $f
        python ~/Scripts/get_seq2.py $f Btarget.txt $f.targetfa
done

for t in *.targetfa; do                                                     #Blast non-anchor species' best hit peptides backa against anchor species (BvA)
        echo $t
        blastp -query $t -db ./$2.prot.db -out $t.BvA.blastout -evalue 1e-6 -max_target_seqs 1 -num_threads $3 -outfmt 6 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 
done

rm Btarget.txt                                                              #Cleanup
rm *targetfa

cat *BvA.blastout > protein.BvA.out                                         #Combine results
rm *BvA.blastout                                                            #Clean up unnecessary files

awk '{print $1, $2}' protein.AvB.out | sort -k 1 > AvB.txt                  #Prepare blastout files of both runs to be in a format that RBH can be retrieved using comm
awk '{print $2, $1}' protein.BvA.out | sort -k 1 > BvA.txt

comm -12 AvB.txt BvA.txt > protein.out                                      #Compare both blastout files keeping only the hits where the best hit in AvB is also the best hit in BvA (RBH)

rm AvB.txt                                                                  #Cleanup
rm BvA.txt

for i in *pep; do                                                           #Run through all the pep files produced before. These are all named after the original transcripts and already contain
        ID=$(basename "$i" .pep)                                            #the corresponding anchor species peptide. Read the peptide from the .pep file and pull out all the hits from the RBH (protein.out)
        echo $ID                                                            #file. Ensure the duplicated anchor species ID is removed. The end result is a list of all peptides across species that
        while read p; do                                                    #are putative one to one orthologs of the original transcript input (after which the file is named). 
                echo $p
                grep "$p" protein.out | awk '{print $2}' | sort -u > "$ID"_pepRBHlist.txt
        done < $i
done

mkdir Fail/                                                                 #Prepare to organise the failed RBH runs
mkdir Fail/ProteinRBH/

for pl in *_pepRBHlist.txt; do                                              #For each peptide RBH list, count the number of peptides present. If the peptides are not present in all species then
        echo $pl                                                            #remove to the Fail directory for other use. Only those that are present in all species will continue in the script.
        ID=$(echo $pl | awk -F "_" -v OFS="_" '{print $1,$2}' )
        COUNT=$(wc -l $pl | awk '{print $1}')
        echo $COUNT
        if [ $COUNT -eq $5 ] ; then
                echo "Blast successfully recovered the correct number of orthologs"
		printf '%s\n' "$ID" >> blast.success.txt
        else
                echo "Blast failed to recover the correct number of orthologs"
                mv $pl Fail/ProteinRBH/
                printf '%s\n' "$ID" >> faillist.txt                          #Keep record of those transcripts that failed to produce RBH peptides across all species in the analysis.
        fi
done

blast=$(ls -lR *_pepRBHlist.txt | wc -l)                                    #Record those that were successfully recovered.

if [ -e faillist.txt ]; then
        mv faillist.txt Proteinfails.txt
	FAILNO=$(wc -l Proteinfails.txt | awk '{print$1}')
        echo "$FAILNO runs failed"
else
        echo "No Failures"
fi

rm *db*                                                                        #Cleanup
rm *.faa

for pl in *_pepRBHlist.txt; do                                                 #Create two protein fasta files per peptide RBH list. One will have the species name followed by transcript ID
        echo $pl                                                               #and the other will be in the standard format. One will be used for alignments and the other to recover transcripts
        ID=$( echo $pl | awk -F "_" -v OFS="_" '{print$1,$2}')                 #from across all species.
        python ~/Scripts/get_seq2.py MasterPep.fa $pl "$ID".fa
        python ~/Scripts/get_seq2.py MasterPepID.fa $pl "$ID".faa
done

mkdir BlastFiles/                                                               #Gather all Blastfiles together.
mv *.out BlastFiles/
mv wanted.fa BlastFiles/
mv protein.fa BlastFiles/
mv translated.txt BlastFiles/

rm Master*                                                                      #Cleanup

mkdir PepRBHlists/                                                              #Gather all peptide RBH lists together. Remove the original .pep files.
mv *_pepRBHlist.txt PepRBHlists/
rm *.pep

cp ~/$4/Transcriptomes/*fasta .                                                 #Get all transcriptomes into the WD to prepare for tblastn.

mv MasterCDS.fasta MasterCDS.fas                                                #Rename MasterCDS file to ensure it isn't included in any later loops. Create another MasterCDS file that consists
cat *.fasta > MasterCDS_ID.fas                                                  #of standard ID format.

# MAKE MULTIPLE SEQUENCE FASTA FILES OF THE TRANSCRIPTS

echo "Generating transcript fasta files via tblastn"

for f in *fasta; do                                                             #Create a blast database of each transcriptome
       echo $f
       makeblastdb -in $f -out $f.db -dbtype nucl -parse_seqids
done

for p in *faa; do                                                               #For each of the protein fastas made above, tblastn against all species transcriptomes in order to retrieve
       ID=$(basename "$p" .faa)                                                 #corresponding transcripts
       echo $ID
       for d in *fasta; do
       		tblastn -query $p -db $d.db -max_target_seqs 1 -evalue 1e-10 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -outfmt 6 -out $ID.$d.blastout -num_threads $3
       done
       cat $ID.*.blastout > $ID.out                                             #For each protein fasta, concatatenate all results into a single blastout file. 
       rm *.blastout                                                            #Remove unnecesssary files.
done

for b in *out; do                                                               #For each blastout file (.out), remove peptides and potentially duplicated (anchor) species' IDs to produce
       ID=$(basename "$b" .out)                                                 #CDS ortholog lists.
       echo $ID
       awk '{print $2}' $b | sort -u > "$ID"_CDSlist.txt
done

mkdir Fail/TblastnRun/                                                          #Prepare for failures

for l in *CDSlist.txt; do                                                       #For each CDS list, count the number of recovered transcripts. If the number is the same as the number of species
       ID=$(echo $l | awk -F "_" -v OFS="_" '{print $1,$2}')                    #included in the run, then the transcript anchor species ID is recorded. Else the out files are moved into the Fail directory.
       COUNT=$(wc -l $l | awk '{print $1}')
       echo $ID $COUNT
       if [ $COUNT -eq $5 ] ; then
		echo "Tblastn recovered all related transcripts"
		printf '%s\n' "$ID" >> tblastn.success.txt
       else
               echo "Tblastn failed to recover correct number of transcripts"
               mv $l Fail/TblastnRun/
               mv $ID.out Fail/TblastnRun/
               printf '%s\n' "$ID" >> faillist.txt
       fi
done

if [ -e faillist.txt ] ; then                                                   #If tblastn runs failed, as is relatively likely with the first, faster method of retrieval, then TCHECK = 1
        mv faillist.txt TBlastNfails.txt
	FAILNO2=$(wc -l TBlastNfails.txt | awk '{print$1}' )
        echo "$FAILNO2 tblastn runs failed"
	declare -i TCHECK=1
else
        echo "No Failures"
	declare -i TCHECK=0
fi


rm *fasta*                                                                      #Cleanup


for c in *CDSlist.txt; do                                                       #For the successful CDS lists, produce multiple sequence transcript fastas per initial anchor transcript input.
        echo $c
        ID=$(echo $c | awk -F "_" -v OFS="_" '{print $1,$2}')
        python ~/Scripts/get_seq2.py MasterCDS.fas $c $ID.fasta
done

rm Master*                                                                      #Cleanup

echo "CDS fasta file generation complete"

#TBlastN second run

if [ $TCHECK -eq 1 ]; then                                                      #This is the second TBlastN check. It will only be activated in the cases where some Tblastn runs failed.
        echo "Recovering TblastN failures"                                      
        mkdir TEMP/                                                             #Produce a temporary directory for fasta files awaiting progression.
        mv *.fa TEMP/
        mv *.faa TEMP/
        mv *.fasta TEMP/
	cp ~/$4/Proteomes/*faa .                                                    #Return the proteomes and transcriptomes to the working directory. Master files here are not necessary
        cp ~/$4/Transcriptomes/*fasta .
        rm Master*
	for f in *faa; do                                                             #For each species, produce a list of just protein IDs.
                species=$(echo $f | awk -F "_" '{print$1}')
                echo $species
                sed -i -e 's/Lalb/lalb/g' $f
                grep ">" $f | awk -F ">" '{print$2}' | awk '{print $1}' | sort -k 1 > $species.protlist
        done
        for f in *fasta; do                                                         #For each species, produce a list of just transript IDs.
                species=$(echo $f | awk -F "_" '{print$1}')
                echo $species
                sed -i -e 's/Lalb/lalb/g' $f
                grep ">" $f | awk -F ">" '{print$2}' | awk '{print $1}' | sort -k 1 > $species.nucllist
        done
        for o in Fail/TblastnRun/*.out; do                                          #For each failed TBlastN run, split the blastout file into transcript(.nucl) and peptides(.prot)
                ID=$(basename $o | awk -F "." -v OFS="." '{print $1,$2}')
                echo $ID
                awk '{print $1}' $o | sort -u | sort -k 1 > $ID.prot
                awk '{print $2}' $o | sort -u | sort -k 1 > $ID.nucl
                for f in *.faa; do                                                  #Loop through each species. Whilst doing so, produce a targettranscript.species.peptide and nucleotide list,
                        species=$(echo $f | awk '{print substr($1,1,4)}' )          #consisting only of those IDs found within that species. For each species, there could be anywhere between 1 and the 
                        comm -12 $species.protlist $ID.prot > $ID.$species.plist    #number of species included in the analysis peptides/transcripts found. 
                        comm -12 $species.nucllist $ID.nucl > $ID.$species.nlist
                        while read p; do                                            #For each of the peptides found in this species, pull out the corresponding blast hits from the original blast file.
                                echo "$p"                                           #Store each hit in a temporary file.
                                grep "$p" $o > $ID.$p.$species.ptmp
                        done < $ID.$species.plist
                        cat *.ptmp | awk '{print $1,$2}' | sort -k 1 > $ID.$species.protout     #Combine peptide hits
                        rm *ptmp                                                                #Remove temp files
                        while read n; do                                            #Repeat the above for transcripts from this species.
                                echo "$n"
                                grep "$n" $o > $ID.$n.$species.ntmp
                        done < $ID.$species.nlist
                        cat *.ntmp | awk '{print $1,$2}' | sort -k 1 > $ID.$species.nuclout
                        rm *ntmp
                        comm -12 $ID.$species.protout $ID.$species.nuclout > $ID.$species.out       #Compare the two files, transcript vs protein. Only those that are present in both (ie the correct
                done                                                                                #protein/transcript combination for that particular species) are kept.
                cat $ID.*.out > $ID.rec_out                                         #Once this has been done for all species, combine the results into one recovery out file.
       		rm *protout
		rm *nuclout
        done
        for r in *rec_out; do                                                       #For each recovery out file, remove peptides and remove duplicates to produce CDSlists.
                echo $r
                ID=$(echo $r | awk -F "." -v OFS="." '{print $1,$2}')
                echo $ID
                awk '{print$2}' $r | sort -u > "$ID"_CDSlist.txt
        done
        mkdir Fail/Recovery/                                                        #Prepare for failures
        for l in *CDSlist.txt; do                                                   #Count through each CDSlist and check the number of recovered transcripts matches that of the species
                ID=$(echo $l | awk -F "_" -v OFS="_" '{print $1,$2}')               #included in the run. Successful recoveries are recorded, failures are recorded and removed into the fail directory.
                COUNT=$(wc -l $l | awk '{print $1}')
                echo $ID $count
                if [ $COUNT -eq $5 ] ; then
                        echo "TblastN Recovery has recaptured the correct number of transcripts"
			            printf '%s\n' "$ID" >> recovery.success.txt
                else
                        echo "Tblastn Recovery has failed to recapture the correct number of transcripts"
                        mv $l Fail/Recovery/
                        mv $ID.rec_out Fail/Recovery/
                        printf '%s\n' "$ID" >> RecoveryFail.txt
                fi
        done
        if [ -e RecoveryFail.txt ] ; then                                           #If there were failures, state how many. Else proclaim complete success.
                FAILNO3=$(wc -l RecoveryFail.txt | awk '{print$1}')
                echo "$FAILNO3 Recovery attempts failed"
        else
                echo "All CDSlists recovered"
        fi
	    cp ~/$4/Transcriptomes/MasterCDS.fasta .                                   #Bring in the Master CDS fasta file.
        mv MasterCDS.fasta MasterCDS.fa
        for c in *CDSlist.txt; do                                                   #Produce fasta files of each CDS list.
                echo $c
                ID=$(echo $c | awk -F "_" -v OFS="_" '{print $1,$2}')
		        echo $ID
                python ~/Scripts/get_seq2.py MasterCDS.fa $c $ID.fasta
		        mv $ID.fasta TEMP/
        done
        rm Master*                                                                  #Cleanup
        echo "Recovered CDS fasta file generation complete"
        rm *nucllist
        rm *protlist
        rm *plist
        rm *nlist
        rm *.out
        rm *prot
        rm *nucl
        rm *fasta
        rm *.faa
        mkdir Recovery/
        mv *rec_out Recovery/
        mv TEMP/*fa .
        mv TEMP/*faa .
	    mv TEMP/*fasta .
        rm -r TEMP/
else
       echo "TBlastN Recovery Step Not Necessary. Continuing"
fi

# CLEANUP

mkdir PepFastaForAlignment/                                                             #More cleanup
mkdir FastaForAlignment/
mkdir CDSlists/
mv *CDSlist.txt CDSlists/
mv *fails.txt Fail/
mv *fail.txt Fail/

rm *faa

for f in *.fasta; do                                                    #Record all the fasta files that made it this far.
        ID=$(basename "$f" .fasta)                                      #Move the fasta file and its corresponding peptide fasta file to their respective folders.
	echo $ID
        echo "$ID is ready for peptide alignment"
        printf '%s\n' "$ID" >> OrthoSuccess.txt
        mv $f FastaForAlignment/
        mv $ID.fa PepFastaForAlignment/
done


CHECK=$(ls -lR *.fa | wc -l )                                       #Check if any protein fasta files have been left behind (ie have no corresponding transcript fasta files)

if [ $CHECK -gt 0 ] ; then                                          #Move any occurences to the fail directory
        mkdir Fail/ProteinSansCDS/
        echo "$CHECK protein fastas did not recover transcript fastas"
        for f in *.fa; do
                ID=$(basename "$f" .fa)
                mv $f Fail/ProteinSansCDS/
                printf '%s\n' "$ID" >> ProteinSansCDSlist.txt    
        done
else
        echo "No protein fastas are alone"
fi

if [ -e ProteinSansCDSlist.txt ] ; then
        mv ProteinSansCDSlist.txt Fail/ProteinSansCDS/
fi

WIN=$(wc -l OrthoSuccess.txt | awk '{print $1}')
echo "$WIN products are ready for alignment"

# PROTEIN ALIGNMENT

echo "Producing Protein Alignments"

FORALN=$(ls -lR PepFastaForAlignment/*fa | wc -l)                       #Note how many target transcripts are now ready for protein products across number of species to be aligned.
echo "$FORALN protein fasta files are preparing for alignment"

bash ~/Scripts/MafftLoop.sh PepFastaForAlignment fa amino               #Use an outside script to loop through the PepFastaForAlignment directory producing alignments (.aln) and dendrograms.

ALN=$(ls -lR *aln | wc -l)                                              #Note how many were aligned.

mkdir Fail/PepAlignment/                                                #Prepare for failures.

if [ $FORALN -eq $ALN ]; then
	echo 'All protein fasta files were successfully aligned'
else
	DIFF=$(($FORALN - $ALN ))                                                          #If there are any failures, move through the protein fasta files and if they lack a corresponding alignment,
	echo "$ALN alignments were produced. $DIFF protein fastas were not aligned"        #move them into the fail directory.
	for f in PepFastaForAlignment/*fa; do
		ID=$(basename "$f" .aln)
		echo $ID
		if [ -e ./$ID.aln ]; then
			echo "$ID peptides were successfully aligned"
		else
			echo "$ID failed to be properly aligned"
			printf "%s\n" "$ID" >> ProtAlignmentFail.txt
			mv $f Fail/PepAlignment/
		fi
	done
fi

mkdir PepAlignment/                                                 #Gather alignments
mv *aln PepAlignment/

for t in PepFastaForAlignment/*tree; do                             #Rename and gather dendrograms
	ID=$(echo $t | awk -F "." '{print$1,$2}' OFS=".")
	mv $t $ID.tree
done

mkdir PepDendrograms/
mv PepFastaForAlignment/*tree PepDendrograms/

# PAL2NAL

echo "Creating codon alignments"

bash ~/Scripts/Pal2NalLoop.sh PepAlignment/ FastaForAlignment/ Pal2Nal/         #Use an outside script to loop through all corresponding transcript fasta files and protein alignments and produce
                                                                                #codon alignments for use in PAML
COD=$(ls -lR Pal2Nal/*paml | wc -l)                                             #Note successful codon alignments.
for p in Pal2Nal/*paml; do
	ID=$(basename "$p" .paml)
	printf "%s\n" "$ID" >> Pal2NalSuccess.txt
done

if [ -e Pal2Nal/Failed/pal2nalfail.txt ]; then                                  #If there are failures, record and move to the fail directory (fasta file and alignments)
	NOCOD=$(ls -lR Pal2Nal/Failed/*paml | wc -l)
	echo "$COD codon alignments successfully created. $NOCOD failed. It is recommended to try again using the Pal2Nal webserver where possible. (-output paml -nogap)"
	mv Pal2Nal/Failed/ Pal2Nal/Pal2NalFail/
	mv Pal2Nal/Pal2NalFail Fail/
    for p in Fail/Pal2NalFail/*paml; do
        ID=$(basename "$p" .paml)
        cp FastaForAlignment/$ID.fasta Fail/Pal2NalFail/
        cp PepAlignment/$ID.aln Fail/Pal2NalFail/
    done
else
	echo "$COD codon alignments successfully created."
fi


# REPORT FROM RUN                                                                           #An overall report of successful runs.
echo "$IN target transcripts input" >> ScriptReport
echo "$blast transcripts were found as protein RBH across $5 species' transcriptomes" >> ScriptReport
echo "$WIN multiple sequence protein fastas were paired with transcript counterparts" >> ScriptReport
echo "$ALN protein alignments were produced" >> ScriptReport
echo "$COD codon alignments were successfully generated and are ready for use in further analyses" >> ScriptReport
printf "\nLauren Mee \nBarribeau Lab \nUniversity of Liverpool" >> ScriptReport

mkdir Success/
mv *uccess.txt Success/

echo "Good luck with PAML"

