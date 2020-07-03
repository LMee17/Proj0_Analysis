* [Show All Code](#)
* [Hide All Code](#)
* [Download Rmd](#)

# Project Zero

# Project Notes

This project is concerned with looking at evolutionary signatures of selection on immune-associated genes, both canonical and candidate. All immune genes (canon/noncanon) are based on orthologs from Apis mellifera (Amel_HAv3.1, NCBI genome).

Two types of analyses have been previously ran using PAML's module, codeml, using gene product codon alignments consisting of 11 bee species, including Apis mellifera. Gene sets are divided into "immunity" and "random", with immunity consisting of "canon" and "non-canon". All alignments consist of 1-1 single-copy orthologs present in all 11 species. Those that did not include all species were dropped.

Evolutionary rate of each gene product was estimated using a codeml run that considers the whole gene product across the given phylogeny as a whole. It is not expected to get dN/dS ratio (omega) values anywhere near 1 (positive selection) using this method as site-specific signals of selection will be likely lost as the ratio is averaged along the length of the product. Instead, this is expected to offer a relative scale of evolutionary rate.

The main analysis consisted of branch-site specific analyses of selection wherein the phylogenetic tree was designated a target, splitting the tree into foreground and background. Codeml then estimates omega value for foreground branches relative to background. For the purposes of this project, each of the designations correspond with the origin or elaboration of social living, to assess if and how these changes in lifestyle affected selection on immune and candidate immune genes in the bees, relative to other genes (the "random" set, which encompass all orthologs that could be made from all proteins in Apis mellifera using the ortholog bash script, OrthoScript [imaginitively named]).

In the null version of this test, omega ratios on the foreground branches are fixed at 1, with the alternative version starting with omega estimates above 1 (a direct test of positive selection). Each model results in a global maximum log likelihoood score which can then be compared using a LRT.

LRT = (lnLAlt - lnLNull)*2

LRT scores can then be compared to a chi-square distribution with 1 degree of freedom to assess significance. LRT scores above 3.84 are considered significant.

P-values will be adjusted using during this analysis to correct for multiple testing using the Benjamini-Hochburg procedure.

# Setting the Verse

In order to annotate the data with class of gene, description of gene (as it is described in Apis mellifera), etc. we first need a "dataverse" of Amel_HAv3.1.

```
source("verseSet.R", local = knitr::knit_global())
```

GeneVerse: created by reading in a file with all genes from Amel HAv3.1 complete with gene product description and transcript ID.

```
head(gene.verse)
```

```
  GeneID   TranscriptID                             Gene_Desc
1    Ant NM_001010975.1                  ADP/ATP translocase
2    Wat NM_001011562.1  worker-enriched antennal transcript
3    Jhe NM_001011563.1            juvenile hormone esterase
4  Mrjp8 NM_001011564.2          major royal jelly protein 8
5    Oa1 NM_001011565.1                  octopamine receptor
6  Kr-h1 NM_001011566.1                    kruppel homolog 1
```

ImmVerse: created by reading in a file with all immune - associated genes used as the input for the orthologue finding program which created the codon alignments used in the codeml branch-site selection analyses, complete with classifications, transcript, peptide and gene IDs.

```
head(imm.verse)
```

```
    TranscriptID    GeneID      ProteinID    Class
1 NM_001011563.1       Jhe NP_001011563.1 NonCanon
2 NM_001011570.1     Vha16 NP_001011570.1    Canon
3 NM_001011574.1 LOC406081 NP_001011574.1 NonCanon
4 NM_001011578.1        Vg NP_001011578.1 NonCanon
5 NM_001011579.1     Mrjp1 NP_001011579.1 NonCanon
6 NM_001011592.1    Tspan6 NP_001011592.1    Canon
```

