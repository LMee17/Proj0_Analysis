#Set the verse so that input codeml result files can be annotated


#Geneverse - created by reading in a file with all genes from Amel HAv3.1 complete with gene 
#product description and transcript ID.
#gene.verse <- read.table("Genome_Misc/Amel_HAv3.1_CDS_GeneDesc.txt", header = T, sep ="\t", quote="")
#head(gene.verse)

#Immverse - created by reading in a file with all immune - associated genes used as the input
#for the orthologue finding program which created the codon alignments used in the codeml 
#branch-site selection analyses, complete with classifications, transcript, peptide and gene 
#IDs.
#imm.verse <- read.table("ImmResources/ImmInput_CodeML_ProjZero_Mar2020.tsv", header =T, sep = "\t")
#head(imm.verse)

#Nov2021
verse <- read.table("Genome_Misc/AlnInfo_ImmFunc_Nov2021.tsv",
                    header = T, sep = "\t")
head(verse)
