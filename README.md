# MAPtools: Command-Line Tools for Mapping-by-Sequencing and QTL-Seq Analysis and Visualization

MAPtools is a set of command-line tools for the analysis and visualization of mapping-by-sequencing (MBS) and QTL-seq data. A typical pipeline involves: (1) aligning the reads of an MBS or QTL-seq experiment using a read mapper, such as Bowtie2 or BWA, (2) processing the alignment files with variant calling software, such as BCFtools or GATK, and (3) analyzing the resulting VCF (Variant Call Format) files with the commands of MAPtools.

## Features

MAPtools v1.00 includes the following commands:

* ``mbs``  performs the analysis of MBS data
* ``qtl``  performs the analysis of QTL-seq data
* ``merge`` re-analyzes the output of `mbs` or `qtl` by binning the allele counts of multiple markers, focusing on haplotypes rather than individual markers
* ``mbsplot``  creates plots using the output of ``mbs`` and `merge`
* ``qtlplot``  creates plots using the output of ``qtl`` and `merge`
* ``annotate`` evaluates the effect of SNPs and indels within a user-specified interval
* ``citation`` gives citation information for MAPtools

## Installation

MAPtools is a Python3 application that can be downloaded from the [MAPtools GitHub page](https://github.com/hcandela/MAPtools), as follows:

```
git clone https://github.com/hcandela/MAPtools.git
```

MAPtools' dependencies can be installed with pip:

```
pip install docopt numpy scipy pandas biopython matplotlib
```

## Contributing

If you would like to contribute to MAPtools, please contact us at hcandela@umh.es.
All contributions are welcome and appreciated!

## License

MAPtools is distributed under the GPL v3.0 License. See LICENSE for more information

## Getting started

### Install dependencies and third-party software:

MAPtools has been tested in pipelines involving the following software:

* **Bowtie2** version **2.4.2**
* **BWA** version **0.7.17-r1188**
* **SAMtools** and **BCFtools** version **1.16**
* **GATK** version **4.0.5.1**

## Tutorial

We provide a step-by-step protocol to reproduce the analyses presented in Figure 1 of the MAPtools article. In this figure, we re-analyze publicly available QTL-seq data from Castillejo *et al*. (2020) and MBS data from Viñegra de la Torre *et al*. (2022).

> Castillejo *et al*. (2020). Allelic variation of MYB10 is the major force controlling natural variation in skin and flesh cholor in strawberry (*Fragaria* spp.) fruit. *Plant Cell* **32:** 3723-3749

> Viñegra de la Torre *et al*. (2022). Flowering repressor AAA+ ATPase 1 is a novel regulator of perennial flowering in *Arabis alpina*. *New Phytologist* **236:** 729-744

### Analysis of QTL-seq data

For this tutorial, you will need to download the raw sequencing data from Castillejo *et al*. (2020), which are available from the [NCBI SRA database](https://www.ncbi.nlm.nih.gov/sra) with accession numbers ERR4463153, ERR4463154, ERR4463155, and ERR4463156. You will also need to download a FASTA file with the [reference genome sequence of *Fragaria vesca*](https://www.rosaceae.org/species/fragaria_vesca/genome_v4.0.a1).

#### 1. Map the reads to the reference genome using your read mapper of choice

1.1. With Bowtie2:

Index the reference genome using bowtie2-build:

```
bowtie2-build Fragaria_vesca_v4.0.a1.fasta.gz INDEX
```

Map the reads to the reference genome using bowtie2. Repeat this step for each available sample:

```
bowtie2 -x INDEX --no-mixed --no-discordant --no-unal -1 ERR4463153_1.fastq.gz -2 ERR4463153_2.fastq.gz  -S ERR4463153.sam
```

1.2. With BWA:

Index the reference genome using the index command of BWA:

```
bwa index -p INDEX Fragaria_vesca_v4.0.a1.fasta.gz 
```

Map the reads to the reference genome using the mem command of BWA. Repeat this step for each available sample:

```
bwa mem INDEX ERR4463153_1.fastq.gz ERR4463153_2.fastq.gz -o ERR4463153.sam
```

#### 2. Convert the SAM files into sorted BAM files using samtools

```
samtools sort -o ERR4463153.bam ERR4463153.sam
```

In this particular experiment, each bulk has been sequenced twice. For this reason, we merge the BAM files for the ERR4463155 and ERR4463156 samples, both of which with white fruits:

```
samtools merge -o white.bam ERR4463155.bam ERR4463156.bam
```

```
samtools index white.bam
```

Similarly, merge the BAM files for the ERR4463153 and ERR4463154 samples, both of which with red fruits:

```
samtools merge -o red.bam ERR4463153.bam ERR4463154.bam
```

```
samtools index red.bam
```

#### 3. Perform the variant calling using BCFtools or GATK

3.1 With BCFtools:

Use the `mpileup` and `call` commands of BCFtools to prepare a single VCF file containing the AD (Allele Depth) field for each bulk:

```
bcftools mpileup -f Fragaria_vesca_v4.0.a1.fasta --annotate FORMAT/AD white.bam red.bam | bcftools call -mv -o variants.vcf
```

You might need to index the genome with ``samtools faidx``

3.2 With GATK:

If you plan to do the variant calling with GATK, note that GATK requires BAM files with the RG field set. You can add this field with the `addreplacerg` command of samtools:

```
samtools addreplacerg -r ID:red -r SM:red -o fixed_red.bam red.bam

samtools addreplacerg -r ID:white -r SM:white -o fixed_white.bam white.bam
```

Use the `HaplotypeCaller` command of GATK to prepare a single VCF file containing the AD (Allele Depth) field for each bulk:

```
gatk HaplotypeCaller --output variants.vcf -I fixed_white.bam -I fixed_red.bam --reference Fragaria_vesca_v4.0.a1.fasta

```

#### 4. Running MAPtools:


4.1 Analyze your data with the ``qtl`` command:

Use the ``qtl`` command to process the ``variants.vcf`` file produced in the previous step. You need to specify which samples are available using the ``--data`` option, indicating them in the same order as in the VCF file. In this experiment, the white and red bulks have been designated as L (low) and H (high), respectively. With the ``-I`` option, the analysis skips the indels and focuses on nucleotide substitutions:

```
./maptools.py qtl --data L,H -i variants.vcf -o qtl.txt -I
```

4.2 Generate plots with the ``plot`` command:

Use the ``plot`` command to automatically generate plots using the output of the previous command, which in this case was saved to the ``qtl.txt`` file:

```
./maptools.py plot -i qtl.txt --captions -m -c Fvb1,Fvb2,Fvb3,Fvb4,Fvb5,Fvb6,Fvb7 -A 800 -a --bonferroni --ci95 -O jpg -o graphics_qtl/
```

In the above example, we used the following options:

``--captions`` generates figure legends
``-m`` produces multi-panel plots
``-c`` selects specific chromosomes to plot
``-A`` specifies how many markers to include in the moving average
``-a`` generates all possible plots
``--bonferroni`` adds Bonferroni threshold to p-value plots
``--ci95`` marks the 95% confidence interval in delta plots
``-o`` specifies the output file format (default is pdf)
``-o`` specifies the output folder

### Analysis of MBS data

For this tutorial, you will need to download the raw sequencing data from Viñegra de la Torre *et al*. (2022), which are available from the [NCBI SRA database](https://www.ncbi.nlm.nih.gov/sra) with accession numbers SRR15564670, SRR11140832, and SRR11140833. Available sequenced samples include a bulk of mutant individuals from an isogenic segregating population (SRR15564670), and two replicates of the resequencing of the non-mutagenized parent involved in the cross (SRR11140832, SRR11140833).

You will also need to download the FASTA file with the [reference genome sequence of *Arabis alpina*](http://www.arabis-alpina.org/data/ArabisAlpina/assemblies/V5.1/), and the [genome annotation file in GFF3 format](http://www.arabis-alpina.org/data/ArabisAlpina/assemblies/V5.1/)

#### 1. Map the reads to the reference genome using your read mapper of choice

1.1. With Bowtie2:

Index the reference genome using bowtie2-build:

```
bowtie2-build Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz INDEX
```

Map the reads to the reference genome using bowtie2. Repeat this step for each available sample:

```
bowtie2 -x INDEX --no-mixed --no-discordant --no-unal -1 SRR15564670_1.fastq.gz -2 SRR15564670_2.fastq.gz  -S SRR15564670.sam
```

1.2. With BWA:

Index the reference genome using the index command of BWA:

```
bwa index -p INDEX Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz 
```

Map the reads to the reference genome using the mem command of BWA. Repeat this step for each available sample:

```
bwa mem INDEX SRR15564670_1.fastq.gz SRR15564670_2.fastq.gz -o SRR15564670.sam
```

#### 2. Convert the SAM files into sorted BAM files using samtools

```
samtools sort -o mutant.bam SRR15564670.sam

samtools sort -o SRR11140832.bam SRR11140832.sam

samtools sort -o SRR11140833.bam SRR11140833.sam
```

Because the parent has been sequenced twice, merge the BAM files for the SRR11140832 and SRR11140833 samples:

```
samtools merge -o parent.bam SRR11140832.bam SRR11140833.bam
```

```
samtools index parent.bam
```

```
samtools index mutant.bam
```

#### 3. Perform the variant calling using BCFtools or GATK

3.1 With BCFtools:

Use the `mpileup` and `call` commands of BCFtools to prepare a single VCF file containing the AD (Allele Depth) field for each bulk:

```
bcftools mpileup -f Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz  --annotate FORMAT/AD mutant.bam parent.bam | bcftools call -mv -o variants.vcf
```

You might need to index the genome with ``samtools faidx``

3.2 With GATK:

If you plan to do the variant calling with GATK, note that GATK requires BAM files with the RG field set. You can add this field with the `addreplacerg` command of samtools:

```
samtools addreplacerg -r ID:mutant -r SM:mutant -o fixed_mutant.bam mutant.bam

samtools addreplacerg -r ID:parent -r SM:parent -o fixed_parent.bam parent.bam
```

Use the `HaplotypeCaller` command of GATK to prepare a single VCF file containing the AD (Allele Depth) field for each bulk:

```
gatk HaplotypeCaller --output variants.vcf -I fixed_mutant.bam -I fixed_parent.bam --reference Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz

```

#### 4. Running MAPtools:

4.1 Analyze your data with the ``mbs`` command:

Use the `mbs` command to process the ``variants.vcf`` file produced in the previous step. You need to specify which samples are available using the ``--data`` option, indicating them in the same order as in the VCF file. In this experiment, a bulk of phenotypically recessive individuals and a wild-type parent have been designated as R and Pd, respectively.

```
./maptools.py mbs --data R,Pd -m R --EMS --parental-filter -o mbs.txt -i variants.vcf
```

In the above example, we used the following options:

``--data`` indicates the available samples
``-m R`` indicates that the mutation is recessive
``--EMS`` focuses on G/C-to-A/T transition mutations
``--parental-filter`` excludes variants also present in the Pd sample
``-o`` specifies the output file

4.2 Generate plots with the ``plot`` command:

Use the ``plot`` command to automatically generate plots using the output of the previous command, which in this case was saved to the ``mbs.txt`` file:

```
./maptools.py plot -i mbs.txt --captions -a -m -c chr1,chr2,chr3,chr4,chr5,chr6,chr6,chr7,chr8 -A 10 --bonferroni -O jpg -o graphics_mbs/
```

In the above example, we used the following options:

``--captions`` generates figure legends
``-a`` generates all possible plots
``-m`` produces multi-panel plots
``-c`` selects specific chromosomes to plot
``-A`` specifies how many markers to include in the moving average
``--bonferroni`` adds Bonferroni threshold to p-value plots
``--ci95`` marks the 95% confidence interval in delta plots
``-O`` specifies the output file format (default is pdf)
``-o`` specifies the output directory

4.3 Annotate the variants with the ``annotate`` command:

Use the ``annotate`` command to predict the effect of the variants detected. It takes the options specified in the previous analysis (mbs or qtl). In the mbs example, we have used the ``-R`` option to focus on an interval located on chromosome 8.

```
./maptools.py annotate -i mbs.txt -g gff_arabis.gff3 -f Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz -R chr8:14250000-16700000 -o annotate_mbs.txt
```

In the above example, we used the following options:

``-i`` input file produced by the maptools mbs or qtl analysis
``-g`` GFF annotation file
``-f`` reference genome file
``-R`` chromosome and interval to analyze
``-o`` specifies the output file

## **Citation**

If you find this software useful, please cite:

César Martínez-Guardiola, Ricardo Parreño, Héctor Candela (2023).
MAPtools: Command-Line Tools for Mapping-by-Sequencing and QTL-Seq Analysis and Visualization. Submitted.
