#!usr/bin/bash

#USAGE bash <script> <MasterIDtable> <isoreducedproteinlist>

#Make a list of the genes involved
for p in *paml; do
        gene=$(echo $p | awk -F "." '{print$1}')
        printf "%s\n" "$gene" >> gene.list
        sort -u gene.list > genereduced.list
done

gcount=$(wc -l genereduced.list | awk '{print $1}')

#Skip genes were there are only one transcript variant
while read g; do
        echo $g
        tcount=$(grep -c "$g" gene.list)
        if [ $tcount == 1 ]; then
                echo "Only one transcript variant found for $g. Skipping."
                if [ -e Output/ ]; then
                        mv $g*paml Output/
                else
                        mkdir Output/
                        mv $g*.paml Output/
                fi
        else
                echo "$tcount transcript variants found for $g. Considering..."
                #First we need to make sure that they are all the same file repeated....
                for i in $g*paml; do
                        echo $i
                        printf "%s\n" "$i" >> $g.transcript.tmp
                done 
                one=$(head -n 1 $g.transcript.tmp)
                grep -v "$one" $g.transcript.tmp > $g.compare.tmp                       #Separate one transcript from the others to allow one vs all comparison
                while read t; do
                        diffc=$(diff $one $t | wc -l)
                        if [ $diffc == 0 ] ; then
                                echo "....."
                        else 
                                echo "Warning. $g alignments are not identical."
                                printf "%s\n" "$g" >> NonIdenticalAlignments.txt
                                if [ -e NonIdentical/ ]; then
                                        mv $g* NonIdentical/
                                else
                                        mkdir NonIdentical/
                                        mv $g* NonIdentical/
                                fi
                        fi
                done < $g.compare.tmp
                awk -F "." '{print $2,$3}' OFS="." $g.transcript.tmp > $g.nucl.tmp
                while read n; do
                        pep=$(grep $n $1 | awk '{print $3}')
                        pepc=$(grep -c $pep $2)
                        if [ $pepc == 1 ]; then
                                echo "$n is the preferred transcript variant"
                                if [ -e Output/ ]; then
                                        mv $g.$n.paml Output/
                                else
                                        mkdir Output/
                                        mv $g.$n.paml Output/
                                fi
                                leftover=$(find -maxdepth 1 -name "$g*.paml" | wc -l)
                                if [ $leftover == 0 ]; then
                                        continue
                                else
                                        if [ -e FilteredOut/ ]; then
                                                mv $g*paml FilteredOut/
                                        else
                                                mkdir FilteredOut/
                                                mv $g*paml FilteredOut/
                                        fi
                                fi
                                break
                        else
                                echo "$n is not the preferred transcript"
                                if [ -e FilteredOut/ ]; then
                                        mv $g.$n.paml FilteredOut/
                                else
                                        mkdir FilteredOut/
                                        mv $g.$n.paml FilteredOut/
                                fi
                        fi
                done < $g.nucl.tmp
                #check that there is still at least one gene file in the Output file
                if [ -e Output/$g*.paml ]; then
                        echo "$g multiple transcript variants have been resolved"
                else
                        echo "$g's preferred transcript variant is not present. Please consider separately."
                        printf "%s\n" "$g" >> Fail.txt
                fi
        fi
done < genereduced.list

rm *tmp
rm *list

fcount=$(ls -lR Output/*paml | wc -l)

if [ $gcount == $fcount ]; then
        echo "All $gcount genes were filtered successfully"
else
        dcount=$(($gcount - $fcount))
        echo "$dcount genes failed to be filtered. Please check Fail.txt"
fi

