# GwasWA Manual

## Contents

- [GwasWA Manual](#gwaswa-manual)
  - [Contents](#contents)
- [Installation](#installation)
  - [Install according to file](#install-according-to-file)
  - [or download the library from conda](#or-download-the-library-from-conda)
- [General parameters](#general-parameters)
- [Module1 WGS data processing](#module1-wgs-data-processing)
    - [1. Download sequencing data, -- step downloadsra](#1-download-sequencing-data----step-downloadsra)
    - [2. sra to fastq,-- step sratofastq](#2-sra-to-fastq---step-sratofastq)
    - [3. fastq quality control, -- step readsqc](#3-fastq-quality-control----step-readsqc)
    - [4. Quality assessment,-step qualityevaluation](#4-quality-assessment-step-qualityevaluation)
    - [5. Establish a reference genome index, -- step downloadref](#5-establish-a-reference-genome-index----step-downloadref)
    - [6. comparison of reference genome, -- step align](#6-comparison-of-reference-genome----step-align)
    - [7, processing bam files, -- step dealbam](#7-processing-bam-files----step-dealbam)
    - [8. Mutation detection, -- step detect](#8-mutation-detection----step-detect)
    - [9、jointgenotype，--step jointgenotype](#9jointgenotype--step-jointgenotype)
    - [10. vcf quality control, -- step vcfqc](#10-vcf-quality-control----step-vcfqc)
- [Module2 GWAS pre-processing](#module2-gwas-pre-processing)
    - [Genotype filling, -- step impute](#genotype-filling----step-impute)
    - [vcf to bfile,-- step transvcf](#vcf-to-bfile---step-transvcf)
    - [gwas quality control, -- step gwasqc](#gwas-quality-control----step-gwasqc)
    - [pca analysis, -- step pca](#pca-analysis----step-pca)
    - [kinship analysis, -- step kinship](#kinship-analysis----step-kinship)
- [Module3 Association analysis](#module3-association-analysis)
    - [Association analysis, -- step association](#association-analysis----step-association)
    - [Screening for significant variation,--step selectsnp](#screening-for-significant-variation--step-selectsnp)
- [Module4 Assessment of variant functional effect](#module4-assessment-of-variant-functional-effect)
  - [Variation impact assessment,--step assess](#variation-impact-assessment--step-assess)
- [Quick Start](#quick-start)
  - [Module 1 WGS data Processing, using E. coli WGS data SRR1770413 as an example](#module-1-wgs-data-processing-using-e-coli-wgs-data-srr1770413-as-an-example)
    - [Download sequencing data, -- step downloadsra](#download-sequencing-data----step-downloadsra)
    - [sra to fastq, -- step sratofastq](#sra-to-fastq----step-sratofastq)
    - [fastq quality control, -- step readsqc](#fastq-quality-control----step-readsqc)
    - [Quality evaluation, --step qualityevaluation](#quality-evaluation---step-qualityevaluation)
    - [Download the reference genome and build its index, --step downloadref](#download-the-reference-genome-and-build-its-index---step-downloadref)
    - [Compare reference genome, -- step align](#compare-reference-genome----step-align)
    - [Process bam files, -- step dealbam](#process-bam-files----step-dealbam)
    - [Variant detection, -- step detect](#variant-detection----step-detect)
    - [jointgenotype，--step jointgenotype](#jointgenotype--step-jointgenotype)
    - [vcf quality control, -- step vcfqc](#vcf-quality-control----step-vcfqc)
  - [Modules 2 and 3](#modules-2-and-3)
    - [vcf to bfile,-- step transvcf](#vcf-to-bfile---step-transvcf-1)
    - [gwas quality control, -- step gwasqc](#gwas-quality-control----step-gwasqc-1)
    - [pca analysis, -- step pca](#pca-analysis----step-pca-1)
    - [kinship analysis, -- step kinship](#kinship-analysis----step-kinship-1)
    - [Association analysis, -- step association](#association-analysis----step-association-1)
    - [Screening for significant variation,-step selectsnp](#screening-for-significant-variation-step-selectsnp)
  - [Module 4 Assessment of variant effect, using human non-coding variant rs11644125 as an example](#module-4-assessment-of-variant-effect-using-human-non-coding-variant-rs11644125-as-an-example)
    - [Variant effect assessment, -- step assess](#variant-effect-assessment----step-assess)

# Installation

## Install according to file

`conda env create -f environment.yml`

`pip install -r requirements.txt`

## or download the library from conda

```bash
conda create -n pipe python=3.9 -y
conda install tensorflow==2.10.0 -y
conda install tensorflow-hub==0.12.0 -y
conda install kipoiseq==0.7.1c -y
conda install matplotlib==3.7.1 -y
conda install beautifulsoup4==4.11.2 -y
conda install seaborn==0.12.2 -y
conda install joblib==1.2.0 -y
conda install scikit-learn==1.0.2 -y
conda install ncbi-datasets-cli==14.15.0 -y
conda install bwa==0.7.17 -y
conda install trim-galore==0.6.10 -y
conda install samtools==1.6 -y
conda install multiqc==1.14 -y
conda install gatk4==4.3.0.0 -y
conda install plink==1.90b6.21 -y
conda install admixture==1.3.0 -y
conda install gcta==1.93.2beta -y
conda install r-base==4.2.2 -y
conda install r-data.table==1.14.8 -y
conda install r-ggplot2==3.4.0 -y
conda install r-reshape2==1.4.4 -y
conda install r-pheatmap==1.0.12 -y
conda install r-qqman==0.1.8 -y
conda install gemma==0.98.3 -y
conda install ensembl-vep==109.3 -y
conda install perl-compress-raw-zlib==2.202 -y
conda install perl-bioperl==1.7.8 -y
conda install htslib==1.9 -y
conda install -c bioconda pyfaidx -y
```

Gwaswa can be accessed globally by adding the gwaswa folder to the environment variable

`export PATH="/path/to/gwaswa:$PATH"`

Add variable settings at the end of the file \~/.bashrc. Next, execute the following command to make the variable settings take effect:

`source ~/.bashrc`


# General parameters

-   \-o, -- output \<path>. Set the directory of the output file, the current directory by default.
-   \-- nosave. Intermediate files are not saved.
-   \-- core \<int>. The number of processes running at the same time.
-   \-- nMem \<str>. Maximum memory footprint.
-   \-- nThrds \<int>. Number of multiple threads.

# Module1 WGS data processing

### 1. Download sequencing data, -- step downloadsra

-   \-- sra \<str>. Downloads based on one or more sra accession, with spaces separating sra accession.
-   \-- sralist \<filename>. Download according to srr\_list.txt. Each line in srr\_list.txt is an sra accession.
-   \-- core \<int>. Number of simultaneous downloads.

`gwaswa --step downloadsra --sra SRR1111111 [SRR2222222 ...]`

`gwaswa --step downloadsra --sralist srr_list.txt`

The sample.sra file downloaded according to the sra accession is stored in the gwaswaOutput/wgs/sra directory. If the download fails, the err\_sra\_log.txt file is generated in the gwaswaOutput/wgs/sra directory to store the failed sra accession.

### 2. sra to fastq,-- step sratofastq

-   \-- sradir \<path>. Provides a directory that holds every sra file that needs to be converted to fastq format.
-   \-- core \<int>. the number of processes that convert core sra files into fastq format at the same time, including compression of fastq format.
-   \-- nThrds \<int>. Number of threads, number of threads used by the process.

`gwaswa --step sratofastq --sradir gwaswaOutput/wgs/sra `

After the sample.sra file in the input directory is converted, the resulting sample.fastq file is stored in the gwaswaOutput/wgs/raw directory in the. gz format.

### 3. fastq quality control, -- step readsqc

-   \-- rawfastqdir \<path>. Provides a directory for each fastq file that needs to be quality controlled.
-   \-- quality \<int>. The default Phred value is 20. The Phred threshold is set and the linker sequence is removed while the low quality bases at the 3 'end are trimmed.
-   \-- phred \<str>. Select the Illumina version to use, optional:
    -   phred33: default. Applies to illumina 1.9: Instruct cutadapt to use an ASCII 33 quality score as the par score.
    -   phred64: Applies to Illumina 1.5: Instruct cutadapt to use the ASCII 64 quality score as the par score.
-   \-- length \<int>. The default value is 20,0 means that this option is not set. set the length threshold. if the length of read is less than this threshold after quality control cleaning or uncoupling, it will be rejected. For double-ended results, if one read in a pair of reads is discarded for this reason, the corresponding other read is also discarded. will not be output to the double ended result file. A value of 0 disables this behavior.
-   \-- stringency \<int>. The default value is 1. How many bases of the linker sequence can be allowed to remain at the end.
-   \-- error \<float>. The default value is 0.1. Maximum allowable error rate.
-   \-- core \<int>. The number of processes that perform fastq quality control at the same time, including compression of fastq format.
-   \-- nThrds \<int>. Number of threads used by each process

`gwaswa --step readsqc --rawfastqdir gwaswaOutput/wgs/raw`

The fastq file after quality control will be stored in the gwaswaOutput/wgs/clean directory in the compressed format of .fq.gz.

### 4. Quality assessment,-step qualityevaluation

-   \-- fastqdir\<path>. Provide a directory of fastq files that require quality assessment.
-   \-- nThrds \<int>. Number of threads used for quality assessment

`gwaswa --step qualityevaluation --fastqdir gwaswaOutput/wgs/clean`

The quality assessment file contains the results of the quality assessment using fastqc and multiqc and is stored in the gwaswaOutput/wgs/qualityEvaluation directory.

### 5. Establish a reference genome index, -- step downloadref

-   \-- accession \<str>. If a local reference genome file is not available, an NCBI Reference sequence can be provided to download the reference genome sequence.

    `gwaswa --step downloadref  --refaccession GCF_000001735.4`
-   \-- taxon \<str>. If a local reference genome file is not available, the NCBI Taxonomy ID or taxonomy name can be provided to download the reference genome sequence.
-   \-- refgenome \<filename>. Use this parameter if you have a local reference genome file.

    `gwaswa --step downloadref --refgenome example/wgs/ref.fa.gz`
-   \--indexalgorithm \<str>。
    -   Is: Fast, but requires a large memory. The database cannot be used when it is larger than 2G.
    -   bwtsw: By default, the sequence needs to be greater than or equal to 10MB to use, which can be used to build larger genomic data.

The reference genome sequence is downloaded and stored in the gwaswaOutput/wgs/ref directory, and the reference genome index file is stored in the same directory as the reference genome.

### 6. comparison of reference genome, -- step align

-   \-- cleanfastqdir \<path>. Provide a directory to store each fastq file after quality control.
-   \-- alignalgorithm\<str>. The algorithm used to align the reference genome.
    -   mem: default. This algorithm is recommended when the reads length is in the range of 70bp-1Mbp.
    -   bwasw: when reads has frequent gap, it is more sensitive than, and this algorithm is recommended. Reads are typically 70bp-1Mbp in length.
    -   backtrack: If the reads length is less than 70bp, this algorithm is recommended. It is recommended that the reads length is less than 100bp. backtrack contains:
        -   aln: aln command to align individual reads to a reference sequence
        -   samse: then use samse or sampe to generate the sam file.
        -   sampe: Use samse or sampe to generate the sam file.
-   \-- refgenome \<filename>. Local reference genome files, use this parameter.
-   \-- core \<int>. Number of fastq files to be compared at the same time.
-   \-- nThrds \<int>. Compare the number of threads used per fastq file

`gwaswa --step align --cleanfastqdir gwaswaOutput/wgs/clean --refgenome gwaswaOutput/wgs/ref/ref.fa `

After genome comparison of the fastq file in the input directory, the sample.bam file is generated in the gwaswaOutput/wgs/align directory.

### 7, processing bam files, -- step dealbam

-   \-- bamdir \<path>. Provide a directory that stores each post-alignment bam file.
-   \-- refgenome \<filename>. Use this parameter if you have a local reference genome file.
-   \-- delPCR. Remove pcr duplicates
-   \-- core \<int>. The number of bam files to be processed at the same time.
-   \-- nThrds \<int>. Number of threads used to process each bam file

`gwaswa --step dealbam --bamdir gwaswaOutput/wgs/align --refgenome gwaswaOutput/wgs/ref/ref.fa`

The bam files in the input directory are sorted, pcr duplicates are removed, indexes are built, and sample\_marked.bam, sample\_marked.bam.bai are generated in the gwaswaOutput/wgs/processed directory.

### 8. Mutation detection, -- step detect

-   \-- processedbamdir \<path>. Provide a directory where each processed bam file is stored.
-   \-- refgenome \<filename>. Use this parameter if you have a local reference genome file.
-   \-- core \<int>. The number of bam files with simultaneous mutation detection.
-   \-- nThrds \<int>. Detects the number of threads used per bam file

`gwaswa --step detect --processedbamdir gwaswaOutput/wgs/processed --refgenome gwaswaOutput/wgs/ref/ref.fa`

After the bam file in the input directory is detected, sample\_g.vcf and its index file are generated in the gwaswaOutput/wgs/gvcf directory.

### 9、jointgenotype，--step jointgenotype

-   \-- gvcfdir \<path>. Provide a directory where each gvcf file is stored.
-   \-- refgenome \<filename>. Use this parameter if you have a local reference genome file.
-   \-- core \<int>. Number of gvcf files split by chromosome at the same time.

`gwaswa --step jointgenotype --gvcfdir gwaswaOutput/wgs/gvcf --refgenome gwaswaOutput/wgs/ref/ref.fa`

In order to reduce the jointgenotype time, first, the\_g.vcf file of each sample in the input directory is divided by chromosome and stored in the gwaswaOutput/wgs/gvcf\_chr directory. Secondly, merge all samples by chromosome, and get chrX\_g.vcf and its index file in gwaswaOutput/wgs/vcf directory. Then, the chrX\_g.vcf is re-compared respectively to obtain the chrX\_vcf file. Finally, the chrX\_vcf files are merged to generate genotype.vcf and its index files, which are stored in the gwaswaOutput/wgs/vcf directory.

### 10. vcf quality control, -- step vcfqc

-   \-- vcfdir \<filename>. Provide the vcf file that holds the variant genotype information.
-   Do hard filtering for SNPs
    -   \-- snpQUAL \<float>. The default value is 30.0. The variation quality value is the value of QUAL in the VCF and is used to measure the reliability of the variation.
    -   \-- snpQD \<float>. SNPQualByDepth, the default value is 2.0. QD is the ratio of the variant quality value (Quality) divided by the depth of coverage (Depth).
    -   \-- snpMQ \<float>. RMSMappingQuality, the default value is 40.0. It is more accurate to describe the degree of dispersion of the quality value of the comparison than the average value.
    -   \-- snpFS \<float>. FisherStrand, the default value is 60.0. The value converted from the p-value of Fisher's test is to describe whether there is obvious positive and negative strand specificity for read containing only variation and read containing only reference sequence bases during sequencing or alignment.
    -   \-- snpSOR \<float>. StrandOddsRatio, the default value is 3.0. SOR was calculated using the symmetric odds ratio test statistical test, corrected for strand specificity.
    -   \-- snpMQRankSum \<float>. MappingQualityRankSumTest, the default value is -12.5.
    -   \-- snpReadPosRankSum \<float>. ReadPosRankSumTest, the default value is -8.0.
-   do hard filtering for indel
    -   \-- indelQUAL \<float>. The default value is 30.0.
    -   \-- indelQD \<float>. SNPQualByDepth, the default value is 2.0.
    -   \-- indelFS \<float>. FisherStrand, the default value is 60.0.
    -   \-- indelSOR \<float>. StrandOddsRatio, the default value is 3.0.
    -   \-- indelMQRankSum \<float>. MappingQualityRankSumTest, the default value is -12.5.
    -   \-- indelReadPosRankSum \<float>. ReadPosRankSumTest, the default value is -8.0.

`gwaswa --step vcfqc --vcfdir gwaswaOutput/wgs/vcf/genotype.vcf --refgenome gwaswaOutput/wgs/ref/ref.fa`

The imported genotype.vcf file undergoes quality control to generate genotype\_filter.vcf and its index file, which is stored in the gwaswaOutput/wgs/vcf directory.

# Module2 GWAS pre-processing

### Genotype filling, -- step impute

-   \-- genotypefile \<filename>. Provide the vcf file that holds the variant genotype information.
-   \-- nMem \<str>. Maximum memory footprint
-   \-- nThrds \<int>. Number of multiple threads used for genotype fill

`gwaswa --step impute --genotypefile gwaswaOutput/wgs/vcf/genotype_filter.vcf.gz --nThrds 10 --nMem 300g`

The input vcf file is filled with genotypes, and the genotype.vcf.gz file is generated in the gwaswaOutput/gwas/transvcf directory.

### vcf to bfile,-- step transvcf

-   \-- genotypefile \<filename>. Provide the vcf file that holds the variant genotype information.
-   \-- phenotypefile \<filename>. Provide the phenotype file. Three columns are required: sample id,family id, and phenotype value, separated by spaces.

`gwaswa --step transvcf --genotypefile gwaswaOutput/gwas/transvcf/genotype.vcf.gz --phenotypefile pheno.txt`

The input vcf file is converted to generate a bfile file and stored in the part2/transvcf directory. Contains. bim,. fam,. bed files. The phenotype file is added to the fam file.

### gwas quality control, -- step gwasqc

-   \-- bfiledir \<path>. Provides a directory where the bfile files are stored.
-   \-- atgc keeps only ATGC
-   \-- snpmiss \<float>. The default value is 0.2. The missing SNPs in the majority of the subjects were excluded. In this step, SNPs of low genotype are deleted.
-   \-- indmiss \<float>. The default value is 0.2. Individuals with high rate of genotype deletion were excluded. In this step, individuals of low genotype are removed. In this step, SNPs of low genotype are deleted.
-   \-- maf \<float>. The default value is 0.05. Minimum allele frequency. SNPs with low MAF are rare and therefore lack the ability to detect SNP phenotypic associations. These SNPs are also more prone to genotyping errors.
-   \-- hwe \<str>. The default value is 1e-6. Exclusion of indicators that deviated from Hardy-Weinberg equilibrium in the control group.
-   \-- hweall \<str>. The default value is 1e-6. Exclude all indicators of sample deviation from Hardy-Weinberg equilibrium.
-   \--indep \<str>。\<window size>\['kb'] \<step size (variant ct)> \<VIF threshold>
-   \--indepPairwise\<str>。\<window size>\['kb'] \<step size (variant ct)> \<r^2 threshold>
-   \--indepPairphase \<str>。\<window size>\['kb'] \<step size (variant ct)> \<r^2 thresh>
-   \-- heterozygosity \<float>. The default value is 3. Exclude individuals with high or low heterozygosity
-   \-- checksex check gender differences.
-   \--rmproblemsex, delete the problem gender
-   \--imputesex, gender filled based on genotype information

`gwaswa --step gwasqc --bfiledir gwaswaOutput/gwas/transvcf`

the bfile file of the input directory is subjected to quality control, and the bfile file after quality control and the intermediate file of quality control are generated and stored in the gwas/qc directory.

### pca analysis, -- step pca

-   \-- cleanbfiledir\<path>. Provides a directory where the bfile files are stored.
-   \-- pcanum \<int>. The default is 6. Number of principal component analyses.
-   \-- groupnum \<int>. By default, the number of groups with the lowest cv error value is selected from 2-20 groups as the number of groups.

`gwaswa --step pca --groupnum 3 --cleanbfiledir gwaswaOutput/gwas/gwasqc`

After group structure analysis and pca analysis, the bfile in the input directory generates pca eigenvalue pca.eigenval and eigenvector pca.eigenvec. The principal component analysis diagram pca.png and group structure diagram admixture.png are drawn and stored in the part2/pca directory.

Principal component analysis figure pca.png

![](image/pca_M2xq0hgdmZ.png)

Group Structure Chart admixture.png

![](image/gt_3_admixture_cFjYSZa-tw.png)

### kinship analysis, -- step kinship

-   \-- cleanbfiledir \<path>. Provides a directory where the bfile files are stored.

`gwaswa --step kinship --cleanbfiledir gwaswaOutput/gwas/gwasqc`

After kinship analysis, the bfile in the input directory generates kinship.txt, draws kinship.png, and stores it in the part2/kinship directory.

kinship.png

![](image/kinship_1Qp28-KldQ.png)

# Module3 Association analysis

### Association analysis, -- step association

-   \-- cleanbfiledir \<path>. Provides a directory where the bfile files are stored.
-   Association analysis model, optional:
    -   \-- lm. generalized linear model.

        `gwaswa --step association --cleanbfiledir gwaswaOutput/gwas/gwasqc --lm`
    -   \-- lmm. mixed linear model.
        -   \-- pcafile \<filename>. Optionally, provide the pca result file as a covariate.
        -   \-- kinshipfile \<filename>. Optionally, provide the kinship result file as a covariate. If it is not provided, it will be automatically generated.
            `gwaswa --step association --cleanbfiledir gwaswaOutput/gwas/gwasqc --lmm --pcafile gwaswaOutput/gwas/pca/pca.eigenvec`

After Association analysis, the bfile in the input directory generates a result.assoc.txt file that stores the information of each variant site, and a variant.vcf file that filters out significantly associated SNPs, draws man.png of Manhattan map and qq.png of qq map, and stores them in the part2/kinship directory.

Manhattan Figure man.png

![](image/man_cFJ0axqfXa.png)

QQ plot qq.png

![](image/qq_JhswH0aLpZ.png)

### Screening for significant variation,--step selectsnp

-   \-- assocfile, association analysis result file. result.assoc.txt
-   \-- pvaluelimit \<str>. Default 1e-5. Filters out snps greater than pvaluelimit.

`gwaswa --step selectsnp --assocfile gwaswaOutput/gwas/association/lm/result.assoc.txt --pvaluelimit 1e-5`

# Module4 Assessment of variant functional effect

## Variation impact assessment,--step assess

\-- species \<str>. Research target species name.

\-- snpfile \<filename>. The association analysis results file variant.vcf is provided as input.

-   SIFT(Sorting Intolerant From Tolerant) is a computational method used to predict the impact of protein variation, which can predict whether a protein variation has a functional impact based on the alignment sequence and protein structure information.

    In SIFT's prediction results, it will be divided into the following four categories according to the degree of influence of variation:
    1.  Tolerated (tolerable): The mutation has no significant effect on the function of the protein, and the mutation site is usually a highly variable region, which is unlikely to be the cause of the disease.
    2.  Tolerated\_Low\_Confidence (low confidence tolerable): The variation has no obvious effect on the function of the protein, but due to the limitation of the aligned sequence and protein structure, the confidence of this prediction result is low and needs further verification.
    3.  Deleterious (deleterious): Variations that have a pronounced effect on the function of a protein may result in changes in the structure or function of the protein that can cause disease or other biological effects.
    4.  Deleterious\_Low\_Confidence (low confidence harmful): The variation has a significant impact on the function of the protein, but due to the limitations of the aligned sequence and protein structure, the confidence of this prediction is low and needs further verification.
        It should be noted that the prediction results of SIFT are based on the alignment sequence and protein structure information, and the prediction results may be affected by the quality of the alignment, the limitation of structural information and other factors, so the prediction results need to be further verified.
-   PolyPhen(Polymorphism Phenotyping) is a computational method used to predict the potential impact of amino acid substitutions on protein structure and function. It provides a score for each replacement, reflecting the probability that the replacement may cause damage to protein function.

    In the PolyPhen scoring system, substitutions are divided into the following three categories:
    1.  Possibly Harmful (probably damaging): The substitution is predicted to have a serious effect on protein function, possibly leading to disease or other biological effects.
    2.  Possibly Harmful (possibly damaging): The substitution prediction may have some effect on protein function, but it is not significant or severe enough to have a minor effect on protein function.
    3.  Benign: The replacement is predicted to have no significant effect on protein function, is usually found in regions of high variation, and is unlikely to cause disease.
        It should be noted that the PolyPhen score is only based on protein structure and function prediction, and its prediction results need to be further verified.

`gwaswa --step vep --snpfile MypipeOutput/gwas/association/lm/variant_lm.vcf`

[^1]: Jiang K, Yang Z, Cui W, et al. An exome-wide association study identifies new susceptibility loci for age of smoking initiation in African-and European-American populations\[J]. Nicotine and Tobacco Research, 2019, 21(6): 707-713.

[^2]: Marees A T, de Kluiver H, Stringer S, et al. A tutorial on conducting genome‐wide association studies: Quality control and statistical analysis\[J]. International journal of methods in psychiatric research, 2018, 27(2): e1608.

# Quick Start

## Module 1 WGS data Processing, using E. coli WGS data SRR1770413 as an example

### Download sequencing data, -- step downloadsra

`gwaswa --step downloadsra --sra SRR1770413 --core 5 --output coli`

The file SRR1770413.sra will be stored in the coli/gwaswaOutput/wgs/sra/path

### sra to fastq, -- step sratofastq

`gwaswa --step sratofastq --sradir coli/gwaswaOutput/wgs/sra --core 5 --nThrds 5 --output coli`

The converted file will be stored in the compressed format of coli/gwaswaOutput/wgs/raw

### fastq quality control, -- step readsqc

`gwaswa --step readsqc --rawfastqdir coli/gwaswaOutput/wgs/raw --core 5 --nThrds 5 --output coli`

The fastq file after quality control will be stored in the coli/gwaswaOutput/wgs/clean directory in the compressed format of .fq.gz.

### Quality evaluation, --step qualityevaluation

`gwaswa --step qualityevaluation --fastqdir coli/gwaswaOutput/wgs/clean --nThrds 5 --output coli`

Contains the results of quality evaluations using fastqc and multiqc, stored in the gwaswaOutput/wgs/qualityEvaluation directory.

### Download the reference genome and build its index, --step downloadref

`gwaswa --step downloadref --accession GCF_000005845.2 --output coli`

The reference genome index file is stored in the same directory as the reference genome.

### Compare reference genome, -- step align

`gwaswa --step align --cleanfastqdir coli/gwaswaOutput/wgs/clean --refgenome coli/gwaswaOutput/wgs/ref/ref.fa --core 5 --output coli`

The .sam file is converted to a .bam file and stored in the coli/gwaswaOutput/wgs/align directory.

### Process bam files, -- step dealbam

`gwaswa --step dealbam --bamdir coli/gwaswaOutput/wgs/align --refgenome gwaswaOutput/wgs/ref/ref.fa --core 5 --nThrds 5 --output coli`

The resulting. bam file is stored in the coli/gwaswaOutput/wgs/processed directory.

### Variant detection, -- step detect

`gwaswa --step detect --processedbamdir coli/gwaswaOutput/wgs/processed --refgenome coli/gwaswaOutput/wgs/ref/ref.fa --core 5 --nThrds 5 --output coli`

After the bam file in the input directory is detected, sample\_g.vcf and its index file are generated in the coli/gwaswaOutput/wgs/gvcf directory.

### jointgenotype，--step jointgenotype

`gwaswa --step jointgenotype --gvcfdir coli/gwaswaOutput/wgs/gvcf --refgenome coli/gwaswaOutput/wgs/ref/ref.fa --core 5 --output coli`

In order to reduce the jointgenotype time, first, the\_g.vcf file of each sample in the input directory is divided by chromosome and stored in the gwaswaOutput/wgs/gvcf\_chr directory. Secondly, merge all samples by chromosome, and get chrX\_g.vcf and its index file in gwaswaOutput/wgs/vcf directory. Then, the chrX\_g.vcf is re-compared respectively to obtain the chrX\_vcf file. Finally, the chrX\_vcf files are merged to generate genotype.vcf and its index files, which are stored in the gwaswaOutput/wgs/vcf directory.

### vcf quality control, -- step vcfqc

`gwaswa --step vcfqc --vcfdir coli/gwaswaOutput/wgs/vcf/genotype.vcf --refgenome coli/gwaswaOutput/wgs/ref/ref.fa --output coli`

The imported genotype.vcf file undergoes quality control to generate genotype\_filter.vcf and its index file, which is stored in the gwaswaOutput/wgs/vcf directory.

## Modules 2 and 3 
Modules 2 and 3 deal with GWAS data and association analysis, taking the experiment of Jiang K et al. as an example[^1], with parameters configured according to the study of Marees A T et al[^2].

### vcf to bfile,-- step transvcf

`gwaswa --step transvcf --genotypefile gwaswa/example/genotype.vcf.gz --phenotypefile gwaswa/example/pheno.txt --output example`

The input vcf file is converted to generate a bfile file and stored in the part2/transvcf directory. Contains. bim,. fam,. bed files. The phenotype file is added to the fam file.

### gwas quality control, -- step gwasqc

`gwaswa --step gwasqc --bfiledir example/gwaswaOutput/gwas/gwasqc --snpmiss 0.2 --indmiss 0.2 --output example`

`gwaswa --step gwasqc --bfiledir example/gwaswaOutput/gwas/gwasqc --snpmiss 0.02 --indmiss 0.02 --output example`

`gwaswa --step gwasqc --bfiledir example/gwaswaOutput/gwas/gwasqc --checksex --output example`

`gwaswa --step gwasqc --bfiledir example/gwaswaOutput/gwas/gwasqc --imputesex --output example`

`gwaswa --step gwasqc --bfiledir example/gwaswaOutput/gwas/gwasqc --maf 0.05 --output example`

`gwaswa --step gwasqc --bfiledir example/gwaswaOutput/gwas/gwasqc --hwe 1e-6 --output example`

`gwaswa --step gwasqc --bfiledir example/gwaswaOutput/gwas/gwasqc --hweall 1e-10 --output example`

`gwaswa --step gwasqc --bfiledir example/gwaswaOutput/gwas/gwasqc --indepPairwise 50 5 0.2 --output example`

the bfile file of the input directory is subjected to quality control, and the bfile file after quality control and the intermediate file of quality control are generated and stored in the gwas/qc directory.

### pca analysis, -- step pca

`gwaswa --step pca --groupnum 3 --cleanbfiledir example/gwaswaOutput/gwas/gwasqc --output example`

The results of population structure and principal component analysis are stored in the gwas/pca directory.

### kinship analysis, -- step kinship

`gwaswa --step kinship --cleanbfiledir example/gwaswaOutput/gwas/gwasqc --output example`

The results of the kinship analysis are stored in the gwas/pca directory.

### Association analysis, -- step association

`gwaswa --step association --cleanbfiledir example/gwaswaOutput/gwas/gwasqc --lm --output example`

The association analysis results are stored in the gwas/association directory.

### Screening for significant variation,-step selectsnp

`gwaswa --step selectsnp --assocfile example/gwaswaOutput/gwas/association/lm/result.assoc.txt --pvaluelimit 1e-5 --output example`

Filter results are stored in the gwas/selectsnp directory

## Module 4 Assessment of variant effect, using human non-coding variant rs11644125 as an example

### Variant effect assessment, -- step assess

`gwaswa --step assess --species homo_sapiens --snpfile example/rs11644125.vcf`

vep and enformer results are stored in the assessment directory