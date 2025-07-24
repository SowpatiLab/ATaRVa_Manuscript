# Steps to generate Haploid bams
This document has the steps to generate haplotype-specific BAM files from an assembly FASTA file, which will be used as input to `aam.py` for generating the truth allele set.

## Download & Index
```bash
# Download Assembly fasta 
$ wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.0.1.fasta.gz

# Index the Assembly
$ samtools faidx hg002v1.0.1.fasta.gz
```

## Chromosome-wise fasta file
```bash
# extract and make fasta file every contig id (both paternal & maternal)
$ samtools faidx hg002v1.0.1.fasta.gz chr1_MATERNAL > chr1_MATERNAL.fa.gz

$ samtools faidx hg002v1.0.1.fasta.gz chr1_PATERNAL > chr1_MATERNAL.fa.gz

# Index every extracted fasta files
$ samtools faidx chr1_MATERNAL.fa.gz
$ samtools faidx chr1_PATERNAL.fa.gz
```

## Align with hg38 using minimap2
```bash
# align every paternal contigs against hg38
$ minimap2 -ax asm5 -t 4 GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz chr1_PATERNAL.fa.gz | samtools view -bS -| samtools sort -@ 4 > chr1_PATERNAL.bam

# align every maternal contigs against hg38
$ minimap2 -ax asm5 -t 4 GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz chr1_MATERNAL.fa.gz | samtools view -bS -| samtools sort -@ 4 > chr1_MATERNAL.bam
```

## Merging of haplotype-wise bams
```bash
# merging paternal bams
samtools merge -@ 6 hg002_PAT.bam chr10_PATERNAL.bam chr11_PATERNAL.bam chr12_PATERNAL.bam chr13_PATERNAL.bam chr14_PATERNAL.bam chr15_PATERNAL.bam chr16_PATERNAL.bam chr17_PATERNAL.bam chr18_PATERNAL.bam chr19_PATERNAL.bam chr1_PATERNAL.bam chr20_PATERNAL.bam chr21_PATERNAL.bam chr22_PATERNAL.bam chr2_PATERNAL.bam chr3_PATERNAL.bam chr4_PATERNAL.bam chr5_PATERNAL.bam chr6_PATERNAL.bam chr7_PATERNAL.bam chr8_PATERNAL.bam chr9_PATERNAL.bam chrY_PATERNAL.bam

# merging maternal bams
samtools merge -@ 6 hg002_MAT.bam chr10_MATERNAL.bam chr11_MATERNAL.bam chr12_MATERNAL.bam chr13_MATERNAL.bam chr14_MATERNAL.bam chr15_MATERNAL.bam chr16_MATERNAL.bam chr17_MATERNAL.bam chr18_MATERNAL.bam chr19_MATERNAL.bam chr1_MATERNAL.bam chr20_MATERNAL.bam chr21_MATERNAL.bam chr22_MATERNAL.bam chr2_MATERNAL.bam chr3_MATERNAL.bam chr4_MATERNAL.bam chr5_MATERNAL.bam chr6_MATERNAL.bam chr7_MATERNAL.bam chr8_MATERNAL.bam chr9_MATERNAL.bam chrX_MATERNAL.bam

# sort the bams
$ samtools sort -@ 6 -o hg002_PAT_sorted.bam hg002_PAT.bam
$ samtools sort -@ 6 -o hg002_MAT_sorted.bam hg002_MAT.bam
```