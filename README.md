# data test
## Small datasets to test IARCbioinfo's pipelines
The data consists of next-generation sequencing data (subsets of whole-exome sequencing with reads ~70bp; Illumina Solexa) from a few individuals from the 1000 genomes project, and reference genome files necessary to run several softwares. The data is subseted to the TP53 gene region of chromosome 17.

## Folder structure

```
data_test/
├── BAM
│   ├── BAM_multiple
│   ├── BAM_realigned_17
│   └── test
├── BED
├── FASTQ
├── REF
│   └── ctat_genome_lib_build_dir_TP53
│       ├── __chkpts
│       ├── _dfam_db_prep_chckpts
│       └── ref_genome.fa.star.idx
└── VCF
```

### BAM folder
- The **BAM** folder contains two BAM files at its root: NA06984_N.bam and NA06984_T.bam with their associated indices (.bai files), that can be used to test somatic variant calling softwares requiring normal/tumor pairs (e.g., workflows IARCbioinfo/mutect-nf and IARCbioinfo/strelka-nf).
- The **BAM/BAM_multiple** folder contains 93 bam files (subsets from different individuals), that can be used to test softwares requiring multiple BAM files (e.g., IARCbioinfo/needlestack)

### BED folder
The **BED** folder contains bed files that can be used by some softwares to specify genomic ranges on which to run the analyses (e.g. variant callers)

### FASTQ folder
The **FASTQ** folder contains raw sequencing files that can be used by alignment softwares (e.g., bwa and STAR using IARCbioinfo/alignment-nf and IARCbioinfo/RNAseq-nf). Although all files are nearly identical, they have different names to allow different kinds of tests: 
- NA06984_N.fastq.gz and NA06984_T.fastq.gz for tumor-normal pairs input
- NA06984_T_1.fastq.gz and NA06984_T_2.fastq.gz for separated read pairs input
- NA06984_T_RG1_1.fastq.gz, NA06984_T_RG1_2.fastq.gz, NA06984_T_RG2_1.fastq.gz, NA06984_T_RG2_2.fastq.gz for demultiplexed inputs

### REF folder
The **REF** folder contains reference files of different sizes (17.fasta and TP53.fasta for the entire chromosome 17 or only the TP53 gene region) and indexes for different softwares. 
- files \*.ht2 are aligner hisat2 indices
- files \*.fasta.0123, .fasta.amb, .fasta.ann, .fasta.bwt, .fasta.bwt.2bit.64, .fasta.bwt.8bit.32, .fasta.fai, .fasta.pac, .fasta.sa are bwa and bwa-mem2 indices
- file 17.dict is a dictionary file, for example used by GATK
- files dbsnp_138.17_7572000-7591000_nochr.vcf.gz (and its index .tbi) and dbsnp_138.17_7572000-7591000.vcf.gz and related files are known germline variant calls from the DBSNP database, used for example in GATK worflows, respectively for references where chromosome names are without the "chr" prefix (chromosome 17 is then just named "17") and with the prefix (chromosome 17 is named "chr17")
- The **REF/ctat_genome_lib_build_dir_TP53** corresponds to a subset of the [Trinity Cancer Transcriptome Analysis Toolkit](https://github.com/NCIP/Trinity_CTAT/wiki) (CTAT) reference genome bundle, that can be used to run STAR and STAR-Fusion (e.g., to test IARCbioinfo/RNAseq-nf and IARCbioinfo/RNAseq-fusion-nf), and notably contains annotations (file ref_annot.gtf) and all STAR indices in subfolder **ref_genome.fa.star.idx**


### VCF folder
The VCF folder contains multi-sample variant calls in the VCF format

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Matthieu Foll*    | follm@iarc.fr | Developer to contact for support |
  | Nicolas Alcala    | alcalan@iarc.fr | Developer |
  | Tiffany Delhomme  | | Developer |
