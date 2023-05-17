# **MAPTOOLS: a Command-Line Tool for Mapping-by-Sequencing and QTL-Seq Analysis and Visualization.**

MAPtools is a set of command-line tools designed and tested to work with the output of bcftools. The typical pipeline involves aligning with bowtie2, processing the alignment with samtools and bcftools, and then using MAPtools to analyze the results.

## **Features**

MAPtools consists of several commands that are used in a similar way to the commands in samtools and bcftools. The current version includes the following commands:

* ``mbs`` : analyzes results of mapping-by-sequencing.
* ``mbsplot`` : automatically generates various types of figures (and their legends) from the results of ``mbs``.
* ``qtl`` : analyzes results of qtl-seq.
* ``qtlplot`` : generates figures from the results of ``qtl``, similar to ``mbsplot``.
* ``merge``: reanalyzes data using pseudomarkers created by combining the results of multiple individual markers (similar to counting haplotype alleles instead of marker alleles). The results of ``merge`` can also be represented with ``mbsplot`` and ``qtlplot``, depending on whether they are of type ``mbs`` or ``qtl``.
* ``annotate``: evaluates the effect of SNPs contained within a specified interval.

## **Installation**

To install MAPtools, follow these steps:

```
git clone https://github.com/hcandela/MAPtools.git
```

## **Usage**

To use MAPtools, run one of the following commands:

* ### **Analyze results of mapping-by-sequencing**

  ``maptools.py mbs [options]``
* ### **Generate figures from results of mapping-by-sequencing analysis**

  ``maptools.py mbsplot [options]``
* ### **Analyze results of qtl-seq**

  ``maptools.py qtl [options]``
* ### **Generate figures from results of qtl-seq analysis**

  ``maptools.py qtlplot [options]``
* ### **Reanalyze data using pseudomarkers created by combining the results of multiple individual markers**

  ``maptools.py merge [options]``
* ### **Evaluate the effect of SNPs contained within a specified interval**

  ``maptools.py annotate [options]``

For more information on each command and its options, run the command with the --help option.

## **Contributing**

If you would like to contribute to MAPtools, please follow these steps:
All contributions are welcome and appreciated!

## **License**

MAPtools is distributed under the GPL v3.0 License. See LICENSE for more information.

## **Overview**

MAPtools is a set of tools to facilitate the analysis of mapping-by-sequencing and QTL-seq studies.

## **Getting started**

### Install dependences and third-party software:

Running MAPtools requires data in VCF/BCF format. The variant calling software compatible with MAPtools is listed here, as well as other useful sotfware:

* Bowtie2: `apt-get install bowtie2`
* BWA: `apt-get install bwa`
* htslib, samtools and bcftools: download the latest version from the webpage: [http://www.htslib.org/download]()
* SRA toolkit is needed to download the example data:

  `wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.5/sratoolkit.3.0.5-ubuntu64.tar.gz`

  Decompress SRA toolkit:

  `tar -xzvf sratoolkit.3.0.5-ubuntu64.tar.gz`
* pigz is recommended for compressing/decompressing data: `apt-get install pigz`

### Running MAPtools: QTL-seq example

For this example, we are using sequencing data analyzed in: Castillejo et al. (2020). *Allelic variation of MYB10 is the major force controlling natural variation in skin and flesh cholor in strawberry* (*Fragaria* spp.) fruit. Plant Cell, 32(12): 3723-3749.

#### 1. Move to maptools directory and download *Fragaria vesca* reference genome version 4.0.1 from [https://www.rosaceae.org/]():

``wget https://www.rosaceae.org/rosaceae_downloads/Fragaria_vesca/Fvesca-genome.v4.0.a1/assembly/Fragaria_vesca_v4.0.a1.fasta.gz``

#### 2. Download SRA reads:

``./sratoolkit.3.0.5-ubuntu64/bin/fastq-dump --split-files ERR4463153 ERR4463154 ERR4463155 ERR4463156``

#### 3. Compress fastq reads:

```
for i in ERR4463153 ERR4463154 ERR4463155 ERR4463156;
do pigz -p {number_of_threads} $i\_1.fastq && pigz -p {number_of_threads} $i\_2.fastq;
done
```

gzip can be used as well, but it is slower

#### 4. Create the bowtie/bwa index. For this example, we are using bowtie2:

``bowtie2-build Fragaria_vesca_v4.0.a1.fasta.gz INDEX``

#### 5. Align the SRA reads to the reference genome for each read:

```
for i in ERR4463153 ERR4463154 ERR4463155 ERR4463156;
do bowtie2 -x INDEX -p {number_of_threads} --no-mixed --no-discordant -1 $i\_1.fastq.gz -2 $i\_2.fastq.gz --no-unal -S sam_qtl_$i\.sam;
done
```

#### 6. Process the SAM files into sorted BAM files:

```
for i in ERR4463153 ERR4463154 ERR4463155 ERR4463156;
do samtools sort -@ {threads} -o bam_qtl_$i\.bam sam_qtl_$i\.sam;
done
```

#### 7. Merge and indexing BAM files of the same pool:

White strawberry pool:

```
samtools merge -@ {number_of_threads} bam_merged_white_pool.bam bam_qtl_ERR4463155.bam bam_qtl_ERR4463156.bam
```

```
samtools index bam_merged_white_pool.bam
```

Red strawberry pool:
``samtools merge -@ {number_of_threads} bam_merged_red_pool.bam bam_qtl_ERR4463153.bam bam_qtl_ERR4463154.bam``

```
samtools index bam_merged_rede_pool.bam
```

#### 8. Variant calling (use bcftools v1.16 or higher preferably)

```
bcftools mpileup -f Fragaria_vesca_v4.0.a1.fasta --annotate FORMAT/AD bam_merged_white_pool.bam bam_merged_red_pool.bam | bcftools call -mv -o variant_calling_qtl.vcf
```

#### 9. Use MAPtools to process the BCF or VCF files produced by BCFtools:

```
./maptools qtl --data L,H -i variant_calling_out.vcf -o qtl_output.txt -I
```

``-I`` option is used to ignore indels for this example

```
./maptools qtlplot -i qtl_output.txt --captions -m -c Fvb1,Fvb2,Fvb3,Fvb4,Fvb5,Fvb6,Fvb7 -A 800 -a --bonferroni --ci95 -o ../graphics_qtl/
```

``--captions`` Generates figure captions
``-m`` Multi-chromosome plots
``-c`` List of selected chromosomes for visualization
``-A`` Number of adjacent markers to calculate moving average lines
``-a`` Generates all possible plots
``--bonferroni`` Show bonferroni test line in p-value plots
``--ci95`` Show the 95% confidence interval in delta plots
``-o`` Output folder. Requires writing permission

#### Running MAPtools: MBS-seq example

For this example, we are using sequencing data analyzed in: Viñegra de la Torre et al. (2022). *Flowering repressor AAA+ ATPase 1 is a novel regulator of perennial flowering in Arabis alpina*. New phytologist, 236: 729-744.

Specifically, we are using the eop002 mutant for the analysis and visualization. This is a MBS of a single sample sequencing data. In the publication, authors found the mutation responsible of the phenotype in chromosome 8, gene Aa_G106560.

#### 1. Download and compress SRA reads:

* Mutant reads:

```
./sratoolkit.3.0.5-ubuntu64/bin/fastq-dump --split-files SRR15564670
```

* Parent reads (2 accession numbers):

```
./sratoolkit.3.0.5-ubuntu64/bin/fastq-dump --split-files SRR11140832 SRR11140833
```

* Compression:

```
for i in SRR15564670 SRR11140832 SRR11140833;
do pigz -p {number_of_threads} $i\_1.fq && pigz -p {number_of_threads} $i\_2.fq;
done
```

#### 2. Download the reference genome adn the GFF annotation

```
wget http://www.arabis-alpina.org/data/ArabisAlpina/assemblies/V5.1/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz
```

```
wget http://www.arabis-alpina.org/data/ArabisAlpina/assemblies/V5.1/A.alpina.modified.5.1.gff3.gz
```

#### 3. Follow steps 4 - 6 of the QTL example:

```
bowtie2-build Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz INDEX
```

```
for i in SRR15564670 SRR11140832 SRR11140833;
do bowtie2 -x INDEX -p {number_of_threads} --no-mixed --no-discordant -1 $i\_1.fq.gz -2 $i\_2.fq.gz --no-unal -S sam_mbs_$i.sam;
done
```

```
for i in SRR15564670 SRR11140832 SRR11140833;
do samtools sort -@ {threads} -o bam_mbs_$i\.bam sam_mbs_$i\.sam && samtools index bam_mbs_$i\.bam;
done
```

* Parent reads need to be merged:

```
samtools merge -@ {number_of_threads} bam_merged_pep1_pool.bam bam_mbs_SRR11140832.bam bam_mbs_SRR11140833.bam
```

* Indexing:

```
samtools index bam_mbs_SRR15564670.bam
```

```
samtools index bam_merged_pep1_pool.bam
```

#### 4. Perform variant calling:

```
bcftools mpileup -f Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz --annotate FORMAT/AD bam_mbs_SRR15564670.bam bam_merged_pep1_pool.bam | bcftools call -mv -o variant_calling_mbs.vcf
```

#### 5. Running MAPtools:

```
./maptools.py mbs --data R,Pd -m R --EMS --parental-filter -o mbs_out.txt -i variant_calling_mbs.vcf
```

``--data`` Defines the type of data. In this case we gave to the mpileup the recesive sample pool (R) and the parental dominant (Pd)
``-m R`` Indicates which pool carries the mutant phenotype
``--EMS`` Filter out SNPs not caused by EMS
``--parental-filter`` Filter out variants present in the parental sequencing.
``-o`` Output file

#### 6. visualizing data MAPtools:

```
./maptools.py mbsplot -i mbs_out.txt -m --captions -c chr1,chr2,chr3,chr4,chr5,chr6,chr6,chr7,chr8 -A 10 -b 20 -a -o ../graphics_mbs/
```

``-m`` Multi-chromosome plots
``--captions`` Generates figure captions
``-c`` List of selected chromosomes for visualization
``-A`` Number of adjacent markers to calculate moving average lines
``-b`` Generates boost plot with an moving average line of the given INT (in this case, 20)
``-a`` Generates all possible plots

#### 7. Analyze MBS data with MAPtools:

As we can see in the plots, the mutation is located in chromosome 8 and approximately between positions 14250000 and 16700000, so we can analyze this region:

```
./maptools.py annotate -d R, Pd -m R --EMS --parental-filter -g gff_arabis.gff3 -f Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz -R chr8:14250000-16700000 -o annotate_mbs.txt
```

## Generic MBS analysis

Use bowtie2 to map the reads from each pool to the reference genome:

```
bowtie2-build genome.fasta GENOMEINDEX
```

```
bowtie2 -x GENOMEINDEX -U dominant_pool.fastq.gz --no-unal -S dominant.sam
```

```
bowtie2 -x GENOMEINDEX -U recessive_pool.fastq.gz --no-unal -S recessive.sam
```

Use the sort command of SAMtools to sort the output files by coordinate and convert them to BAM format:

```
samtools sort -o dominant.bam dominant.sam
```

```
samtools sort -o recessive.bam recessive.sam
```

Process the resulting BAM files with the mpileup and call commands of BCFtools:

```
bcftools mpileup -f genome.fasta --annotate FORMAT/AD dominant.bam recessive.bam | bcftools call -mv -V indels -o output.vcf
```

Use MAPtools to process the BCF or VCF files produced by BCFtools:

```
cat output.vcf | maptools mbs -d D,R -r R -m D -o mbs_results.txt
```

```
maptools mbsplot -i mbs_results.txt -c 1,2,3 --captions -A 200 -m --all --bonferroni
```

```
maptools merge -i mbs_results.txt -w 2 -o reprocessed_results.txt
```

```
cat output.vcf | maptools annotate -g genome.gff3 -f reference.fasta -d D,R -r D -m R -o annotation.txt -R 1:1-10000000
```

## **Citation**

César Martínez-Guardiola, Ricardo Parreño, Héctor Candela (2023).
MAPtools: a Command-Line Tool for Mapping-by-Sequencing and QTL-Seq Analysis and Visualization. Submitted.
