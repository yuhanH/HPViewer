# HPViewer: sensitive and specific genotyping of human papillomavirus in metagenomic DNA

##  Description
HPViewer is a tool for <span style="color:blue">__genotyping and quantification of HPV from metagenomic or human genomic shotgun sequencing data.__</span>  We designed it to improve performance by masking nonspecific sequences from reference genomes and directly identifying HPV short DNA reads. It contains two HPV databases with different masking strategies, repeat-mask and homology-mask and one homology distance matrix to choose between those two databases.

__If you use the HPViewer software, please cite our manuscript:__

Hao Y, Yang L, Galvao Neto A, Amin M R, Kelly D, Brown S M, Branski R C, Pei Z. (2017). HPViewer: Sensitive and specific genotyping of human papillomavirus in metagenomic DNA. bioRxiv, 208926. (accepted by Bioinformatics)


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
a) input files (__-U__ or __-1 -2__): fastq files (or fastq.gz), unpaired (-U unpaired.fastq) or R1,R2 paired (-1 R1.fastq -2 R2.fastq)

b) output file name (__-o__)

###  Optional
a) database mask type (__-m__):  hybrid-mask(default), repeat-mask, homology-mask.

If you set -m, it should be in front of reads input (-m repeat-mask -1 R1.fastq -2 R2.fastq). Repeat-mask is a more sensitive mode; and homology-mask is suggested when some types of HPV are present in large abundance which may lead to false positive of other types of HPV.

b) number of threaded used in bowtie2 alignment (__-p__)

c) minimal coverage threshold to determine HPV present (__-c__), default is 150 bp (1.5 x average length of your reads).


##  Results
a) output_HPV_summary.txt has three coloumns with types of HPV present, number of reads per kilobase (RPK) for the matching HPV, and number of reads of the matching HPV.

<p align="center"><img src="http://mathurl.com/yaps6bcj.png" height="30" /></p>



b) alignment results after bowtie2: output.sam, output.bam



##  Basic Usage (demo)
```{bash eval=FALSE}
python HPViewer.py -U test_unpaired.fastq -o TEST
```

```{bash}
more TEST/TEST_HPV_profile.txt
```

##  Work Flow
<p align="center"><img src="https://github.com/yuhanH/HPViewer/blob/master/workflow.png" height="512" /></p>



