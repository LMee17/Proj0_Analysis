#!usr/bin/bash

#USAGE bash <script> <MasterIDtable>

for p in *paml; do
        trans=$(echo $p | awk -F "." '{print$1,$2}' OFS=".")
        echo $trans
        gene=$(grep "$trans" $1 | awk '{print $1}')
        echo $gene
        suffix=$(echo $p | awk -F "." '{print $NF}')
        echo "$suffix"
        mv $p "$gene.$trans.$suffix"
done
