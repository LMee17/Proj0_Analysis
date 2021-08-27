#!usr/bin/bash


#USAGE bash BranchSiteCodeml.sh <run name> <directory of codon alignments> <tree file> <M0 ctl file> <null ctl file> <alt ctl file>

#REQUIREMENTS
#All alignments should end in .paml and be prefixed with gene or protein name followed by a "." (ie Dcr.Aln.paml)
#Tree file should be unrooted with specified branch target(s) included either with #1 or $1
#Directory containing alignments should be within the current working directory


check=$(echo $2 | grep -c "/")
if [ $check -eq 1 ]; then
	WD=$(echo $2 | sed 's/.$//')
	echo $WD
else
	WD=$2
	echo $WD
fi

printf "%s\t%s\t%s\t%s\n" "Gene" "lnL_Null" "lnL_Alt" "LRT" >> $1.lnLResults.txt

tot=$(ls -lR $WD/*paml | wc -l)
echo "$tot alignments given as input to analysis" >> $1.report

now=$(date '+%c')
echo "Analysis began at $now" >> $1.report

sed -i "s/treefile = .*/treefile = $3/" $4

prog=1
for p in $WD/*paml; do
	gene=$(echo $p | awk -F "." '{print $1}' | awk -F "/" '{print$2}')
	file=$(echo $p | awk -F "/" '{print$NF}')
	tik=$(date '+%c')
	echo "Gene $prog/$tot: $gene" >> $1.report
	echo "Began: $tik" >> $1.report
	echo "$prog/$tot: $gene"
	#Run M0
	echo "Running M0 model"
	sed -i "s/seqfile = .*/seqfile = $WD\/$file/" $4
	sed -i "s/outfile = .*/outfile = $gene.$1.M0.out/" $4
	sed -i "s/noisy = .*/noisy = 0 /" $4
	codeml $4
	#Extract tree and other parameters
	#Extract tree
	grep ");" $gene.$1.M0.out | tail -n 1 > $gene.$1.M0.tree
	#();!!
#	sed -i "s/(((Aflo: [0-9]\.[0-9]*, Amel: [0-9]\.[0-9]*): [0-9]\.[0-9]*, ((Bter: [0-9]\.[0-9]*, Bimp: [0-9]\.[0-9]*): [0-9]\.[0-9]*, Mqua: [0-9]\.[0-9]*): [0-9]\.[0-9]*): [0-9]\.[0-9]*, Emex: [0-9]\.[0-9]*)/&\$1/g" $gene.$1.M0.tree		#Corb origin of sociality, longterm shift in selection
#	sed -i "s/(((Aflo: [0-9]\.[0-9]*, Amel: [0-9]\.[0-9]*): [0-9]\.[0-9]*, ((Bter: [0-9]\.[0-9]*, Bimp: [0-9]\.[0-9]*): [0-9]\.[0-9]*, Mqua: [0-9]\.[0-9]*): [0-9]\.[0-9]*): [0-9]\.[0-9]*, Emex: [0-9]\.[0-9]*)/&\#1/g" $gene.$1.M0.tree		#Corbiculate origin of sociality, episodic selection
#	sed -i "s/(((Aflo: [0-9]\.[0-9]*, Amel: [0-9]\.[0-9]*)/&\$1/g" $gene.$1.M0.tree		#Apis elaboration of sociality
#	sed -i "s/, Mqua/&\$1/g" $gene.$1.M0.tree		#Meliponi elaboration of sociality
#	sed -i "s/((Lalb/& \$1/g" $gene.$1.M0.tree		#Halictid origin of sociality
#	sed -i "s/(Ccal/& \$1/g" $gene.$1.M0.tree		#Xylo origin of sociality
	#dN and dS
	dN=$(grep "tree length for dN:" $gene.$1.M0.out | awk '{print $NF}')
	dS=$(grep "tree length for dS:" $gene.$1.M0.out | awk '{print $NF}')
	#kappa
	kappa=$(grep "kappa" $gene.$1.M0.out | awk -F "=" '{print $2}' | xargs)
	#omega
	omega=$(grep "omega" $gene.$1.M0.out | awk -F "=" '{print $2}' | xargs)
	echo "M0 run complete" >> $1.report
	echo "Kappa = $kappa" >> $1.report
	echo "dN = $dN, dS = $dS" >> $1.report
	echo "Omega = $omega" >> $1.report
	#Cleanup
	mkdir $gene.$1/
	mv 2NG.dN $gene.M0.2NG.dN
	mv 2NG.dS $gene.M0.2NG.dS
	mv 2NG.t $gene.M0.2NG.t
	mv 4fold.nuc $gene.M0.4fold.nuc
	mv lnf $gene.M0.lnf
	mv rst $gene.M0.rst
	mv rst1 $gene.M0.rst1
	mv rub $gene.M0.rub
	mv $gene.M0* $gene.$1/
	#Run Null Model
	echo "Running Null model"
	sed -i "s/seqfile = .*/seqfile = $WD\/$file/" $5
	sed -i "s/outfile = .*/outfile = $gene.$1.Null.out/" $5
	sed -i "s/treefile = .*/treefile = $gene.$1.M0.tree/" $5
	sed -i "s/ kappa = .*/ kappa = $kappa /" $5
	sed -i "s/noisy = .*/noisy = 0 "/ $5
	codeml $5
	#Extract Null model lnl
	lnLNull=$(grep "lnL" $gene.$1.Null.out | awk -F ":" '{print $4}' | awk '{print $1}')
	nullNP=$(grep "lnL" $gene.$1.Null.out | awk -F ":" '{print $3}' | grep -o '[0-9]\+')
	echo "Null model run complete. lnL = $lnLNull, np = $nullNP" >> $1.report
	#Cleanup
	mv 2NG.dN $gene.Null.2NG.dN
	mv 2NG.dS $gene.Null.2NG.dS
	mv 2NG.t $gene.Null.2NG.t
	mv lnf $gene.Null.lnf
	mv rst $gene.Null.rst
	mv rst1 $gene.Null.rst1
	mv rub $gene.Null.rub
	mv $gene.Null* $gene.$1/
        #Run Alt Model
	echo "Running Alt models"
	sed -i "s/seqfile = .*/seqfile = $WD\/$file/" $6
	sed -i "s/treefile = .*/treefile = $gene.$1.M0.tree/" $6
	sed -i "s/ kappa = .*/ kappa = $kappa /" $6
	sed -i "s/noisy =.*/noisy = 0 /" $6
        count=1
        #Iterate through omega values until convergence reached
        for o in $(seq 1.2 0.4 3.2); do
                sed -i "s/ omega = .*/ omega = $o /" $6
                sed -i "s/outfile = .*/outfile = $gene.$1.Alt$count.out/" $6
                echo "Running iteration $count, omega = $o" >> $1.report
                codeml $6
                lnLAlt=$(grep "lnL" $gene.$1.Alt$count.out | awk -F ":" '{print $4}' | awk '{print $1}')
                altNP=$(grep "lnL" $gene.$1.Alt$count.out | awk -F ":" '{print $3}' | grep -o '[0-9]\+')
                echo "Alt model run $count complete. lnL = $lnLAlt, np = $altNP" >> $1.report
                #CleanUp
                mv 2NG.dN $gene.Alt$count.2NG.dN
                mv 2NG.dS $gene.Alt$count.2NG.dS
                mv 2NG.t $gene.Alt$count.2NG.t
                mv lnf $gene.Alt$count.lnf
                mv rst $gene.Alt$count.rst
                mv rst1 $gene.Alt$count.rst1
                mv rub $gene.Alt$count.rub
                mv $gene.Alt$count* $gene.$1/
                #Convergence check
                printf "%s\t%s\t%s\n" "$lnLNull" "$lnLAlt" "$count" >> $gene.tmp
                LRT=$(sed "${count}q;d" $gene.tmp | awk '{print ($2-$1)*2; }')
		concheck=$(echo "$LRT" | grep -c "^-")
                if [ $concheck -eq 0 ]; then
                        echo "Convergence appears to have been reached." >> $1.report
                        incheck=0
                        break
                fi
                if [ $o == "3.2" ]; then
                	echo "Max omega iterations used and convergence has continued to fail." >> $1.report
                        incheck=1
		else
			echo "Reiterating ..." >> $1.report
                fi
		let count=count+1
        done
        echo "Checking..." >> $1.report
        if [ $incheck -eq 1 ]; then
                echo "Attempting alt run using in.codeml parameter file." >> $1.report
                grep "lnL" -A 2 $gene.$1.M0.out | sed '1d' | awk '{$1=$1;print}' > in.codeml
		sed -i "s/outfile = .*/outfile = $gene.$1.AltI.out/" $6
		sed -i "s/ omega = .*/ omega = 1.2/" $6
		sed -i "s/ kappa = .*/ kappa = 2/" $6
                codeml $6
                lnLAlt=$(grep "lnL" $gene.$1.AltI.out | awk -F ":" '{print $4}' | awk '{print $1}')
                altNP=$(grep "lnL" $gene.$1.AltI.out | awk -F ":" '{print $3}' | grep -o '[0-9]\+')
                echo "Alt model using in.codeml input complete. lnL = $lnLAlt, np = $altNP" >> $1.report
                printf "%s\t%s\t%s\n" "$lnLNull" "$lnLAlt" "i" >> $gene.tmp
                LRT=$(grep "i" $gene.tmp | awk '{print ($2-$1)*2; }')
                concheck=$(echo "$LRT" | grep -c "^-")
                #CleanUp
                mv in.codeml $gene.inputparameters.txt
                mv $gene.inputparameters.txt $gene.$1/
                mv 2NG.dN $gene.AltI.2NG.dN
                mv 2NG.dS $gene.AltI.2NG.dS
                mv 2NG.t $gene.AltI.2NG.t
                mv lnf $gene.AltI.lnf
                mv rst $gene.AltI.rst
                mv rst1 $gene.AltI.rst1
                mv rub $gene.AltI.rub
                mv $gene.AltI* $gene.$1/
                if [ $concheck -eq 0 ]; then
                        echo "Convergence appears to have been reached. Moving on to the next gene." >> $1.report
                else
                        echo "This gene continues to experience problems." >> $1.report
                        echo "$gene" >> $1.Fail.txt
                        besthit=$(sort -nr -k 2 $gene.tmp | head -n 1 | awk '{print $3}')
                        echo "The best attempt was run $besthit" >> $1.report
                        echo "Logging and moving on" >> $1.report
                        lnLAlt=$(sort -nr -k 2 $gene.tmp | head -n 1 | awk '{print $2}')
                        LRT=$(sort -nr -k 2 $gene.tmp | head -n 1 | awk '{print ($2-$1)*2; }')
                fi
        fi
        #Record
        echo "Recording results..." >> $1.report
        printf "%s\t%s\t%s\t%s\n" "$gene" "$lnLNull" "$lnLAlt" "$LRT" >> $1.lnLResults.txt
        #Final Cleanup
        rm *tmp
        mv *out $gene.$1/
 	mv *tree $gene.$1/
        tok=$(date '+%c')
        echo "Completed: $tok" >> $1.report
        let prog=prog+1
done 

grep -P "\t\t" $1.lnLResults.txt | awk '{print $1}' >> $1.Fail.txt 

mkdir $1/
mv *.$1/ $1/

