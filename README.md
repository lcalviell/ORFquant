# ORFquant
An R package for Splice-aware quantification of translation using Ribo-seq data


*ORFquant* is an R package that aims at detecting and quantifiying ORF translation on complex transcriptomes using Ribo-seq data.
This package uses syntax and functions present in Bioconductor packages like *GenomicFeatures*, *rtracklayer* or *BSgenome*. 
*ORFquant* aims at quantifying translation at the single ORF level taking into account the presence of multiple transcripts expressed by each gene.
To do so, the *ORFquant* pipeline consists of transcript filtering, *de-novo* ORF finding, ORF quantification and ORF annotation.
A variety of annotation methods, both in transcript and genomic space, is performed for each ORF, to yield a more complete picture of alternative splice sites usage, uORF translation, translation on NMD candidates and more.

More details can be found in our manuscript:

### Quantification of translation uncovers the functions of the alternative transcriptome ###

*Lorenzo Calviello^, Antje Hirsekorn, Uwe Ohler^*

**biorXiv (2019)**, doi: https://doi.org/10.1101/608794

https://www.biorxiv.org/content/10.1101/608794v2

We recommend users to have a look at the vignette: https://htmlpreview.github.io/?https://github.com/lcalviell/ORFquant/blob/master/ORFquant_vignette.html, or our manual (*ORFquant_manual.pdf*).


To install *ORFquant*:

```
library("devtools")
install_github(repo = "lcalviell/ORFquant")

library("ORFquant")

```

Three steps are required to use *ORFquant* on your data:
```
?prepare_annotation_files
```
parses a *.gtf* and a *.2bit* file. (this need to be done once per each annotation-genome combination, a .2bit file can be obtained from a fasta file using the *faToTwoBit* software from UCSC: https://genome.ucsc.edu/goldenpath/help/twoBit.html - http://hgdownload.soe.ucsc.edu/admin/exe/ )


```
?prepare_for_ORFquant
```
or (**recommended**) the *Ribo-seQC* package (https://github.com/lcalviell/Ribo-seQC) can create input files for *ORFquant* using a Ribo-seq .bam file.


```
?run_ORFquant
```

is the master function used to perform the entire analysis workflow, for single genes or (**recommended**) entire transcriptomes.
Please check the vignette for an example workflow.


For any question, please email:

calviello.l.bio@gmail.com or uwe.ohler@mdc-berlin.de


Enjoy!


