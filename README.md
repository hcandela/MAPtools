# **MAPTOOLS: a Command-Line Tool for Mapping-by-Sequencing and QTL-Seq Analysis and Visualization.**

MAPtools is a set of command-line tools designed and tested to work with the output of bcftools. The typical pipeline involves aligning with bowtie2, processing the alignment with samtools and bcftools, and then using MAPtools to analyze the results.

## **Features**
MAPtools consists of several commands that are used in a similar way to the commands in samtools and bcftools. The current version includes the following commands:

*   ``mbs`` : analyzes results of mapping-by-sequencing.

*   ``mbsplot`` : automatically generates various types of figures (and their legends) from the results of ``mbs``.

*   ``qtl`` : analyzes results of qtl-seq.
*   ``qtlplot`` : generates figures from the results of ``qtl``, similar to ``mbsplot``.
*   ``merge``: reanalyzes data using pseudomarkers created by combining the results of multiple individual markers (similar to counting haplotype alleles instead of marker alleles). The results of ``merge`` can also be represented with ``mbsplot`` and ``qtlplot``, depending on whether they are of type ``mbs`` or ``qtl``.
*   ``annotate``: evaluates the effect of SNPs contained within a specified interval.

Installation
To install MAPtools, follow these steps:

Usage
To use MAPtools, run one of the following commands:

## Analyze results of mapping-by-sequencing
MAPtools mbs [options] <input_file>

## Generate figures from results of mapping-by-sequencing analysis
MAPtools mbsplot [options] <input_file>

## Analyze results of qtl-seq
MAPtools qtl [options] <input_file>

## Generate figures from results of qtl-seq analysis
MAPtools qtlplot [options] <input_file>

## Reanalyze data using pseudomarkers created by combining the results of multiple individual markers
MAPtools merge [options] <input_file>

# Evaluate the effect of SNPs contained within a specified interval
MAPtools annotate [options] <input_file>
For more information on each command and its options, run the command with the --help option.

Contributing
If you would like to contribute to MAPtools, please follow these steps:
All contributions are welcome and appreciated!

License
MAPtools is distributed under the MIT License. See LICENSE for more information.

Overview
MAPtools is a set of tools to facilitate the analysis of mapping-by-sequencing and QTL-seq studies.

Getting started

Use bowtie2 to map the reads from each pool to the reference genome:

bowtie2-build genome.fasta GENOMEINDEX
bowtie2 -x GENOMEINDEX -U dominant_pool.fastq.gz --no-unal -S dominant.sam
bowtie2 -x GENOMEINDEX -U recessive_pool.fastq.gz --no-unal -S recessive.sam

Use the sort command of SAMtools to sort the output files by coordinate and convert them to BAM format:

samtools sort -o dominant.bam dominant.sam
samtools sort -o recessive.bam recessive.sam

Process the resulting BAM files with the mpileup and call commands of BCFtools:

bcftools mpileup -f genome.fasta --annotate FORMAT/AD dominant.bam recessive.bam | bcftools call -mv -V indels -o output.vcf

Use MAPtools to process the BCF or VCF files produced by BCFtools:

`` cat output.vcf | maptools mbs -d D,R -r R -m D -o mbs_results.txt ``

maptools merge -i mbs_results.txt -w 2 -o reprocessed_results.txt

cat output.vcf | maptools annotate -g genome.gff3 -f genoma.fasta -d D,R -r D -m R -o annotation.txt -R 1:1-10000000


Citation

César Martínez-Guardiola, Ricardo Parreño, Héctor Candela (2023).
MAPtools: a Command-Line Tool for Mapping-by-Sequencing and QTL-Seq Analysis and Visualization.



