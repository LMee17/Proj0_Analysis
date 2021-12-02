#!usr/bin/bash

#USAGE bash CodonCheck.sh < PSG list > < Parsed Results >

#Look at codons individually

#Get Codons per PSG list

for l in $1*txt; do
        run=$(echo $l | awk -F "/" '{print$NF}' | awk -F "." '{print$1}')
        mkdir $run.CodonCheck/
        echo $run
        while read g; do
                printf "%s\t%s\t%s\n" "Codon" "2H:CMD" "3H:CMD" > $run.CodonCheck/$g.CodonCheck.tsv
                awk '$2 > 5 {print ;}' $2$g.2hit.Results.tsv | awk '{print$1}' | sed '1d' | sort -k1 > $g.2hit.tmp
                awk '$2 > 5 {print ;}' $2$g.3hit.Results.tsv | awk '{print$1}' | sed '1d' | sort -k1 > $g.3hit.tmp
                while read c; do
                        codon=$(echo $c | awk '{print$1}')
                        check=$(grep -c -E "^$codon$" $g.2hit.tmp)
                        check3=$(grep -c -E "^$codon$" $g.3hit.tmp)
                        if [ $check -eq 1 ]; then
                                cmd2="TRUE"
                        else
                                cmd2="FALSE"
                        fi
                        if [ $check3 -eq 1 ]; then
                                cmd3="TRUE"
                        else
                                cmd3="FALSE"
                        fi
                        printf "%s\t%s\t%s\n" "$codon" "$cmd2" "$cmd3" >> $run.CodonCheck/$g.CodonCheck.tsv
                done < $run.PSCodons/$g.PSCodons.txt
        done < $l
        printf "%s\t%s\t%s\t%s\n" "Gene" "NoSelectedCodons" "No_2hitCMD" "No_3hitCMD" > $run.CodonCheck.tsv
        for t in $run.CodonCheck/*tsv; do
                gene=$(echo $t | awk -F "/" '{print$NF}' | awk -F "." '{print$1}')
                echo $gene
                nocod=$(sed '1d' $t | wc -l | awk '{print$1}')
                no2cmd=$(sed '1d' $t | awk '{print $2}' | grep -c "TRUE")
                no3cmd=$(sed '1d' $t | awk '{print$3}' | grep -c "TRUE")
        printf "%s\t%s\t%s\t%s\n" "$gene" "$nocod" "$no2cmd" "$no3cmd" >> $run.CodonCheck.tsv
        done
        rm *tmp
done
