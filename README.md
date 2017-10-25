---
title: "Overview of HPViewer"
output:
  rmarkdown::github_document:
    toc: yes
  pdf_document:
    toc: yes
linkcolor: green
---


##  Description
HPViewer is a tool for <span style="color:blue">__genotyping and quantification of HPV from metagenomic or human genomic shotgun sequencing data.__</span>  We designed it to improve performance by masking nonspecific sequences from reference genomes and directly identifying HPV short DNA reads. It contains two HPV databases with different masking strategies, repeat-mask and homology-mask and one homology distance matrix to choose between those two databases.

__If you use the HPViewer software, please cite our manuscript:__

Yuhan Hao, Liying Yang, Antonio Galvao Neto, Milan R. Amin, Dervla Kelly, Stuart M. Brown, Ryan C. Branski, Zhiheng Pei.<span style="color:blue"> __"HPViewer: Sensitive and specific genotyping of human papillomavirus in metagenomic DNA"__</span>  (Submitted)



##  Installation
```{bash eval=FALSE}
$ git clone https://github.com/yuhanH/HPViewer.git
```

###  Pre-requisites

Python (2.7+)

Python packages (sys, getopt, subprocess)

Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

SAMtools: http://www.htslib.org/

Bedtools: http://bedtools.readthedocs.io/en/latest/

##  Parameters
###  Required
a) input files (__-U__ or __-1  -2 __): fastq files (or fastq.gz), unpaired (-U unpaired.fastq) or R1,R2 paired (-1 R1.fastq -2 R2.fastq)

b) output file name (__-o__)

###  Optional
a) database mask type (__-m__):  hybrid-mask(default), repeat-mask, homology-mask. repeat-mask is a more sensitive mode; and homology-mask is suggested when some types of HPV are present in large abundance which may lead to false positive of other types of HPV.
b) number of threaded used in bowtie2 alignment (__-p__)
c) minimal coverage threshold to determine HPV present (__-cov__), default is 150 bp (1.5 x average length of your reads).


##  Results
a) output_HPV_summary.txt has three coloumns with types of HPV present, number of reads per kilobase (RPK) for the matching HPV, and number of reads of the matching HPV.

$$ average\:coverage = \frac{RPK * average\:reads\:length}{1000}\ $$

b) alignment results after bowtie2: output.sam, output.bam



##  Basic Usage (demo)
```{bash eval=FALSE}
python HPViewer2.py -U test_unpaired.fastq -o TEST
```

```{bash}
more TEST/TEST_HPV_profile.txt
```

##  Wrok Flow (demo)
<p align="center"><img src="https://github.com/yuhanH/HPViewer/blob/master/workflow.png" height="256" /></p>



