# Project Zero: Assessing a changing immune system across 11 bee genomes.

## 1: Selection Analyses

### Step 1: Genome Preprocessing 

Genome information for species used in project
- _Apis florea_: Aflo1.0, Kapheim _et al_ (2015)
- _Apis mellifera_: AmelHAv3.1, Wallberg _et al_ (2019)
- _Bombus impatiens_: Bimp1.0, Sadd _et al_ (2015)
- _Bombus terrestris_: Bter1.0, Sadd _et al_ (2015)
- _Ceratina calcarata_: ASM165200v1, Rehan _et al_ (2016)
- _Dufourea novaeangliae_: Dnov1.0, Kapheim _et al_ (2015)
- _Eufriesea mexicana_: Emex1.0, Kapheim _et al_ (2015)
- _Habropoda laboriosa_: Hlab1.0, Kapheim _et al_ (2015)
- _Lasioglossum albipes_: Lalb2.0, Kocher _et al_ (2013)
- _Megachile rotundata_: Mrot1.0, Kapheim _et al_ (2015)
- _Melipona quadrifasciata_: Mqua1.0, Kapheim _et al_ (2015)

All mentioned scripts can be found in GenomeWork/Scripts/.

#### 1.1 Resolving Isoforms

In order to reduce the number of isoforms in each genome, I wrote a script `IsoResolve.sh` which resolves isoforms using, amongst other things, the program `IsoSel` (Phillipon _et al_ 2017, http://doua.prabi.fr/software/isosel). 

To begin with, I parsed the corresponding .gff files of each species using the script `GFFread.sh`. It requires the gff file and proteome as input.

```sh
bash GFFRead.sh <gff file> <proteome>
```

This takes in a GFF and outputs (1) a `MasterIDtable.tsv`, which contains all genes and their respective transcripts and proteins (if present), (2) an `isoforms_locus_tag.txt` file, which is necessary input for `IsoSel` and (3) an `IsoformCount.txt` file, another input file for the next step.

`IsoResolve.sh` requires `IsoSel` and the following as input per species:
- 1: a gene list 
- 2: the proteome
- 3: MasterIDtable.tsv
- 4: IsoformCount.txt
- 5: isoforms_locus_tag.txt
- 6: Output name to be prefixed in all output files

and the python script `get_seq2.py`.

Command (1 line)

```sh
bash IsoResolve.sh AmelGenes.txt \\ 
Amel_GCF_003254395.2_Amel_HAv3.1_protein.faa MasterIDtable.tsv \\
isoforms_locus_tag.txt IsoformCount.txt Amel
```

This script runs by running through a series of logical statements in order to resolve isoforms. First it looks to see if there are annotated isoforms (ie a NP over an XP) to store as the preferred isoform per gene. Should this not be the case and there are multiple isoforms it then recruits IsoSel to attempt to resolve the issue. If IsoSel fails, it then takes the longest isoform instead. If there are multiple annotated isoforms it repeats these steps using just the NP orthlogs. It cannot resolve by transcript, instead relying on the proteome information in order to work. 

As output `IsoResolve.sh` produces an isoform-reduced proteome (in this case, `Amel_Isoform_filtered.faa`), a table with genes and preferred protein ortholog per gene (`Amel.Filteredfaa.table`) and a script run report with other miscellaneous information, including how each gene was resolved.

The isoform-reduced _A.mellifera_ proteome was used in the next data processing step. 

#### 1.2 Determining Orthology

Required input:
- Transcriptomes and proteomes of all species involved (including the isoform reduced _A. mellifera_ proteome)
- List of transcripts of interest

List of transcripts of interest were either a list of canonical immune transcripts, non-canonical immune transcripts or "background" genes. These lists were made using literature (Evans _et al_ 2006, Geer _et al_ 2009, Alaux _et al_ 2011, Richard _et al_ 2012, Doublet _et al_ 2017, Zdobnov _et al_ 2017) and can be found using ImmResources/. 

The script `Orthoscript.sh` is the master script that determines orthology. It works in a number of steps that I would make much simpler if I was writing it for the first time again after learning more about coding!

`OrthoScript.sh` requires the programs `PAL2NAL` (Suyama _et al_ 2006, http://www.bork.embl.de/pal2nal/), `Blast+` (Camacho _et al_ 2009, https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), `MAFFT` (Katoh _et al_ 2017, https://mafft.cbrc.jp/alignment/software/) and the accompanying scripts `get_seq2.py`, `Pal2NalLoop.sh` and `MAFFTLoop.sh`.

The input it requires is:
- 1: A list of transripts. These are transcripts of interest (see above)
- 2: The designated "anchor" species. By this I mean the species which will be used as the primary template against which the other species will be blasted. In this case the anchor species is _Apis mellifera_ as (1) it is the most well annotated species of the builds in this study and (2) the canon / candidate immune transcripts come from experiments that were undergone using this species. It is designated using the prefix of the transcriptome/proteome (see script notes)
- 3: Number of threads for the `BLAST+` steps.
- 4: The working directory. It is expected that within this directory there is a directory named Proteomes/ with the proteomes contained within, and the same for the transcriptomes. Each file has a prefix with the species in the manner laid out in the script's preamble and suffix dependent on whether the fasta file contains protein or nucleotide sequences.
- 5: Number of species included in the run. In this case it is 11.

To run (i.e.)

```sh
bash OrthoScript.sh CanonCDS.txt Amel 20 ~/ProjZero/ 11
```

`OrthoScript.sh` translates the list of transcripts into proteins using `BLAST+` and the isoform-reduced _A. mellifera_ proteome. Orthologs in the other species are identified using reciprocal best-hit `BLAST`, keeping only those where a one-to-one ortholog is found among all species. When all 11 species have a representative protein, multiple sequence protein alignments (MSAs) are generated using `MAFFT`. Each protein ortholog is then reverse translated (using `tBLASTn`) into nucleotide sequences. These MSA/fasta pairs are then fed into `PAL2NAL` which produces codon alignments with gaps removed. 

The end result of `OrthoScript.sh` are codon alignments ready for analysis in `PAML` (Yang 1997, 2007, http://abacus.gene.ucl.ac.uk/software/paml.html). 

#### 1.3 Finishing Steps

The produced alignments are named after the starting _A. mellifera_ transcript (an architect from how I originally envisioned the script to work which would be the first thing I'd remove now) and so I wrote up two scripts to clean up (1) by adding the gene to the file name and (2) ensuring that any multiple isoforms that may have escaped through the net are resolved. (As the transcriptomes were not isoform reduced there were a number of instances were multiple transcripts mapped back to one protein).

The first of these scripts is `GeneName.sh`. It should be run in the directory containing the codon alignments and requires the _A. mellifera_ MasterIDtable.tsv file as input.

The next script to run is then `GeneFilter.sh` which also requires the MasterIDtable.tsv file as input. Further to this, `GeneFilter.sh` requires a list of all the proteins IDs present in the isoform reduced proteome. At the end, all alignments which passed the requirements are output in the folder FilteredFastas/, whereas any that failed are put aside and listed in a fail report for manual checking.

### Step 2: Selection Analyses Using `PAML`

All non-scripting files metioned are found in GenomeWork/PAML_Files/.

#### 2.1 Positive Selection Analyses

To examine whether distinct origins of sociality and social elaboration alter rate of evolution on immunologically relevant genes, we used branch-site likelihood methods for detecting positive selection in `codeml` (Yang and Nielson, 2002, Zhang _et al_ 2005). To run through genes in batches I wrote `BranchSiteCodeml_V2.sh`.

`BranchSiteCodeml_V2.sh` requires `PAML` and the following input
- 1: The name of the run. I.e test or "Eusocial"
- 2: A directory containing the codon alignments that are to be ran through `PAML`. Ideally this should be within the directory where you are running the script. All alignments need to be suffixed ".paml"
- 3: The tree file. In my case this was `anthophila_unrooted_Jan20.newick`. The tree is, perhaps obviously, unrooted and does not contain branch lengths (`PAML` estimates branch lengths per gene analysed).
- 4: M0 ctl file. I used `H0_codeml.ctl`
- 5: The Null model ctl file, `branchsiteNull.ctl`
- 6: The alternative model ctl file, `branchsiteAlt.ctl`

For example: (one line)

```sh
bash BranchSiteCodeml.sh Run1 CodonAlignments/ \\
anthophila_unrooted_Jan20.newick H0_codeml.ctl \\
branchsiteNull.ctl branchsiteAlt.ctl
```

`BranchSiteCodeml.sh` runs through a series of steps per each gene's codon alignment. First, it runs the M0 model (see ctl file for details), which estimates the branch lengths for the gene tree with the fixed topology provided, and the dN, dS, kappa and omega parameters. It then edits the `branchsiteNull.ctl` file, adding in estimated kappa and other details that allow the analysis to move forward (alignment location etc). Once ran, the -lnL is extracted from the output file and stored. The alternative model begins in the same way, except the intial oemga value is edited. For the first run omega = 1.2, the alternative model is allowed to run and -lnL is extracted. If there is a large negative difference between the -lnL of the alternative and null model (i.e. the alternative model is significantly less likely than the null, resulting in a negative LRT score), then the script runs `codeml` again, iterating different values of intial omega in steps of 0.4 upt to a maximum of 3.2. If this still fails, then the script uses the M0 output to make an in.codeml file to dictate all parameters. Should there still be a negative LRT score the most likely iteration is recorded. 

The resulting files from the script run above include the `Run1.report`, which includes all the steps taken per gene during the whole analysis, `Run1.Fail.txt` which keeps track of genes where `codeml` failed altogether, and `Run1.lnLResults.txt`, which consists of each gene, the null log likelihood and alt likelihood and the resultant LRT score.

#### 2.2 Evolutionary Rate

To assess overal evolutionary rate for each gene we used `codeml` to estimate gene-wide _dN/dS_ ratios. For each run, the script `mktable_2.sh` is used to create an output file with all the estimated kappa values and branch-length trees produced by `codeml` per gene in the M0 runs of the above analyses. This output, `Run1.M0.readout` (or whatever the run is named) is used as the input for the next step.

The script `Omega_codeml.sh` requires three inputs:
- 1: the path to the directory of codon alignments (suffixed *paml)
- 2: the ctl file `omega.ctl`
- 3: an M0.readout file such as `Run1.M0.readout`

The script takes the estimated kappa value and phylogeny per gene from the .readout file and adds them into `omega.ctl` file. `codeml` is ran (model = 0, NSsites = 0, ncatG = 1) and the estimated _dN/dS_ ratio is extracted from the results file.  All _dN/dS_ ratios for all genes included in the alignment directory are recorded in the output file `dNdSratio.tsv`. 

### Step 3: FitMultiModel
The branch-site test of positive selection (BST) has been found to be vulnerable to producing false positives when there are instances of multinucleotide mutation events (MNM, Venkat _et al_ 2018). Sites where more than one nucleotide has changed simultaneously (codons of multiple differences, CMDs) violate the model we have used which considers each change as an individual event. CMDs can possibly still have become fixated due to positive selection, else processes of neutral mutation, and the BST cannot distinguish between the two. We decided to highlight any sites identified as likely CMDs as less certain indications of positive selection. 

In order to do this, we enlisted `FitMultiModel` (Lucaci _et al_ 2021) from `HyPhy` to assess our alignments and predicted gene trees and flag any codons likely (those with evidence ratios above 5) to be CMDs using the script `RunNParseFMMBatch.sh`.

`RunNParseFMMBatch.sh` requires the following input
- 1: A list of genes identified as under positive selection. We only ran this program on genes already assessed by paml to be under positive selection.
- 2: A directory containing the same alignments fed into `PAML`. As such they are still assumed to be suffixed ".paml"
- 3: A text file listing the species included in the alignment. Names need to match those used in the alignment.
- 4: A directory containing all the trees predicted by `PAML` per PSG. 

After the script is run there will be two results directories. The first, `ResultFiles/`, contains all the json files produced by `FMM`. The second, `ParsedResults/`, contains .tsv files parsed from the original json output. Per gene there are two tables: a 2Hit and 3Hit results.tsv with codon sites and estimated evidence ratios (ERs). These consider the likelihood of CMDs occurring due to two (2H) or up to three (3H) simultaenous nucleotide changes having occurred. Those with an ER above 5 are considered likely to be CMDs. 

For the next step, the codons `PAML` has designated under selection are needed per PSG per branch test in which it was found under selection. We considered sites that had a Bayes Emperical Bayes posterior probability above 0.9 as likely to be under selection. An easy way of doing this is to cycle through results folders from `PAML` and extracting them thus: 
```sh
sed -n '/^Bayes Empirical Bayes/,/^\The grid/{p;/\n/q}' $p*out | \\
sed '1,2d' | tac | sed '1,3d' | tac | sed 's/\*//g' | \\
awk '$3>0.9 {print;}' > $p.PSCodons.txt
```

where `$p` here referred to a gene under selection from a PSG text list. 

We then ran the script `CodonCheck.sh`. Input consists of two directories: the first containing all the PSG lists for each individual branch test run, the second being the `ParsedResults/` directory produced above. This script results in a results file per branch site test (ie Run1 from above) where each positively selected codon is checked against each gene's likely CMDs, `Run1.CodonCheck.tsv`. A directory, `Run1.CodonCheck/` is generated (per branch test), where all results files are kept per gene. 

### Step 4: Data Analysis

#### 3.1 The Analyses

Each of these steps are explained in the `Proj0.rmd` file.

#### 3.2 Visualisation

All R scripts used in the above .rmd files are found in Scripts/ directory. All plots found in the Plots/ directory were made using the R scripts `MakePlots.R` and `GOPlotPrep.R`.

## 2: Gene Family Analyses

Files used in this analysis are found in GenomeWork/CAFE_Files/.

Paralogs were gathered into gene families based on protein sequence using `fastOrtho` (Davis _et al_ 2020) using the isoform reduced proteomes created as directed above resulting in the output `FastOrtho_countsTable.tsv`. `CAFE5` (Mendes _et al_ 2020), was then was used to assess the statistical significance of gene family change across the phylogeny in light of the sociality of each branch. The tree provided to `CAFE5` with branch lengths was `Tree4_11Bees.newick` and the tree that designated the social lifestyle of each branch was `Tree4_wSociality.newick`. An example of a `CAFE5` run command is below. `-c` and `-o` refer to number of threads and output, respectively.

```sh
cafe5 -in FastOrtho_countsTable.tsv -t Tree4_11Bees.newick \\
-y Tree5_wSociality.newick -c 20 -o Out_Aug21
```

The python notebook scripts, `getCafeTables.ipynb` and `summarizeCafe_3levelsAugust2021.ipynb` were used to assess results and produce results tables and figures.

Please contact Dr Seth Barribeau who underwent this portion of the analysis for furher direction.

## Program Versions

`BLAST` blast 2.7.1, build Oct 18 2017 \
`CAFE5` 5.0.0 \
`FastOrtho` PATRIC3/FastOrtho \
`Linux` 16.04.5 LTS (Xenial Xerus) \
`HyPhy` 2.5.33 \
`MAFFT` 7.407 \
`OrthoFinder` v2.5.4 \
`Pal2Nal` v14 \
`PAML` 4.9 \
`R` 4.1.2 (2021-11-01) -- "Bird Hippie" 

## References
Alaux, C., Dantec, C., Parrinello, H. & Le Conte, Y. Nutrigenomics in honey bees: Digital gene expression analysis of pollen’s nutritive effects on healthy and varroa-parasitized bees. _BMC Genomics_ 12, (2011). \
Camacho C, _et al_. BLAST+: architecture and applications. _BMC Bioinformatics_. (2009) Dec 15;10:421. doi: 10.1186/1471-2105-10-421. PMID: 20003500; PMCID: PMC2803857. \
Davis, J. J. _et al_. The PATRIC Bioinformatics Resource Center: expanding data and analysis capabilities. _Nucleic Acids Research_. Volume 48, Issue D1, 08 January 2020, Pages D606–D612, https://doi.org/10.1093/nar/gkz943 \
Doublet, V. _et al_. Unity in defence: honeybee workers exhibit conserved molecular responses to diverse pathogens. _BMC Genomics_ 18, 207 (2017). \
Emms DM, Kelly S. OrthoFinder: phylogenetic orthology inference for comparative genomics. _Genome Biol_. 20(1):238, (2019) doi: 10.1186/s13059-019-1832-y. PMID: 31727128; PMCID: PMC6857279. \
Evans, J. D. _et al_. Immune pathways and defence mechanisms in honey bees _Apis mellifera_. _Insect Mol Biol_. 15, 645–656 (2006). \
Geer, L. Y. _et al_. The NCBI BioSystems database. _Nucleic Acids Res_. 38, 492–496 (2009). \
Kapheim, K. M. _et al_. Genomic signatures of evolutionary transitions from solitary to group living. _Science_. 348, 1139–1143 (2015). \
Katoh, K., Rozewicki, J. & Yamada, K. D. MAFFT online service: multiple sequence alignment, interactive sequence choice and visualization. _Brief. Bioinform_. 1–7 (2017). doi:10.1093/bib/bbx108 \
Kocher, S. D. _et al_. The draft genome of a socially polymorphic halictid. _Genome Biol_. 14, R142 (2013). \
Lucaci, A. G., Wisotsky, S. R., Shank, S. D., Weaver, S., and Kosakovsky Pond,  S. L. Extra base hits: widespread empirical support for instantaneous multiple-nucleotide changes. _Plos one_, 16(*3*):e0248337, 2021. \
Mendes FK, Vanderpool D, Fulton B, Hahn MW. CAFE 5 models variation in evolutionary rates among gene families. _Bioinformatics_. (2020) doi: 10.1093/bioinformatics/btaa1022. Epub ahead of print. PMID: 33325502. \
Philippon, H., Souvane, A., Brochier-Armanet, C. & Perrière, G. IsoSel: Protein isoform selector for phylogenetic reconstructions. _PLoS One_. 12, 1–13 (2017). \
Rehan, S. M., Glastad, K. M., Lawson, S. P. & Hunt, B. G. The Genome and Methylome of a Subsocial Small Carpenter Bee, _Ceratina calcarata_. _Genome Biol. Evol_. 8, 1401–1410 (2016). \
Richard, F.-J., Holt, H. L. & Grozinger, C. M. Effects of immunostimulation on social behavior, chemical communication and genome-wide gene expression in honey bee workers (_Apis mellifera_). _BMC Genomics_ 13, 558 (2012). \
Sadd, B. M. _et al_. The genomes of two key bumblebee species with primitive eusocial organization. _Genome Biol_. 16, 76 (2015). \
Suyama, M., Torrents, D., Bork, P. & Delbru, M. PAL2NAL : robust conversion of protein sequence alignments into the corresponding codon alignments. _Nucleic Acids Res_. 34, 609–612 (2006). \
Venkat, A., Hahn, M. W. and J. W. Thornton. Multinucleotide mutations cause false inferences of lineage-specific positive selection. _Nature ecology
& evolution_, 2(8):1280–1288, 2018. \
Wallberg, A. _et al_. A hybrid de novo genome assembly of the honeybee, _Apis mellifera_, with chromosome-length scaffolds. _BMC Genomics_. 20, 1–19 (2019). \
Yang, Z. PAML: a program package for phylogenetic analysis by maximum likelihood. _Comput. Appl. Biosci_. 13, 555–6 (1997). \
Yang, Z. PAML 4: Phylogenetic analysis by maximum likelihood. _Mol. Biol. Evol_. 24, 1586–1591 (2007). \
Zdobnov, E. M. Tegenfeldt, F. Kuznetsov, D. Waterhouse, R. M. Simão, F. A. Ioannidis, P. Seppey, M. Loetscher, A. Kriventseva, E. V. OrthoDB v9.1: cataloging evolutionary and functional annotations for animal, fungal, plant, archaeal, bacterial and viral orthologs. _Nucleic Acids Research_. 45(D1) D744 - D749 (2017). https://doi.org/10.1093/nar/gkw1119 \
Zhang, J., Nielsen, R. & Yang, Z. Evaluation of an Improved Branch-Site Likelihood Method for Detecting Positive Selection at the Molecular Level. _Mol. Biol. Evol_. 22, 2472–2479 (2005).












