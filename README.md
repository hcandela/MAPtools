# MAPtools: Command-Line Tools for Mapping-by-Sequencing and QTL-Seq Analysis and Visualization

MAPtools is a set of command-line tools for the analysis and visualization of mapping-by-sequencing (MBS) and QTL-seq data. A typical pipeline involves: (1) aligning the reads of an MBS or QTL-seq experiment using a read mapper, such as Bowtie2 or BWA, (2) processing the alignment files with variant calling software, such as BCFtools or GATK, and (3) analyzing the resulting VCF (Variant Call Format) files with the commands of MAPtools.

## Contributing

If you would like to contribute to MAPtools, please contact us at hcandela@umh.es.
All contributions are welcome and appreciated!

## License

MAPtools is distributed under the GPL v3.0 License. See LICENSE for more information

## Features

MAPtools v1.00 includes the following commands:

* ``mbs``  performs the analysis of MBS data
* ``qtl``  performs the analysis of QTL-seq data
* ``merge`` re-analyzes the output of `mbs` or `qtl` by binning the allele counts of multiple markers, focusing on haplotypes rather than individual markers
* ``plot``  creates plots using the output of ``mbs``, ``qtl`` or `merge`
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

## Getting started

### Install dependencies and third-party software:

MAPtools has been tested in pipelines involving the following software:

* **Bowtie2** version **2.4.2**
* **BWA** version **0.7.17-r1188**
* **SAMtools** and **BCFtools** version **1.16**
* **GATK** version **4.0.5.1**

## Workflow overview

A typical MBS or QTL-seq analysis comprises several steps, including: 

#### 1. Quality trimming of the reads
As a first step, you may want to process the FASTQ files containing the reads with programs such as [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) or [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) to trim adaptor sequences and discard low-quality bases.

#### 2. Indexing the reference genome and read mapping
The next step is to align the reads of each sample to the reference genome sequence. We have tested MAPtools in pipelines involving the [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [BWA](https://github.com/lh3/bwa) read mappers. However, it should also work with any other read mapper, provided that it outputs files in standard SAM or BAM formats. SAM files must be converted into sorted BAM files using [Samtools](http://www.htslib.org/download/). 

In the example, a reference genome is first indexed using `bowtie2-build`:
```
bowtie2-build reference_genome.fa INDEX
```
The reads are then aligned to the reference genome using `bowtie2`, and the output is directly converted to sorted BAM format using `samtools sort`:
```
bowtie2 -x INDEX --no-mixed --no-discordant --no-unal -1 reads_forward.fq.gz -2 reads_reverse.fq.gz | samtools sort -o aligned_reads.bam -
```

#### 3. Variant calling
The alignment files (in BAM format) are then processed altogether using variant calling software to identify variant sites. We have tested MAPtools in pipelines involving [GATK](https://gatk.broadinstitute.org/hc/en-us) and [BCFtools](http://www.htslib.org/download/), but it should also work with any other variant calling software that can output uncompressed VCF format files containing an AD (allelic depth) field for each sample.

In our example, three different BAM files from a mapping-by-sequencing experiment are processed using a pipeline involving `bcftools mpileup` and `bcftools call` to generate a single VCF file (`variants.vcf`):
```
bcftools mpileup -f reference_genome.fasta --annotate FORMAT/AD dominant_bulk.bam recesive_bulk.bam parent.bam | bcftools call -m -v -o variants.vcf
```
With the `-v` flag present, the `variants.vcf` file will only save data for the variant positions. Non-variant positions are also discarded by MAPtools, but larger VCF files increase the processing time.

As an alternative, you may also choose to store the output as a binary BCF file (rather than VCF) by setting the `-O b` option:
```
bcftools mpileup -f reference_genome.fasta --annotate FORMAT/AD dominant_bulk.bam recesive_bulk.bam parent.bam | bcftools call -m -v `-O b` -o variants.bcf
```

#### 4. Running MAPtools
Once your `variants.vcf` file is ready, you can go ahead and run MAPtools. For additional details, please also check the Tutorial section below.

Following on with our example, we run the `mbs` command of MAPtools, as follows:
```
./maptools.py mbs --data D,R,Pd -m R -o mbs.txt -i variants.vcf
```
Here, we use ``--data`` to specify the samples that are available for the analysis, **in the same order** as they were processed while making the VCF file. In our case, these samples are designated as D (the bulk of phenotypically dominant individuals), R (the bulk of phenotypically recessive individuals) and Pd (resequencing of the parent of the mapping population showing the dominant phenotype). We also use ``-m R`` to specify that the mutation is recessive.

If you previously chose to store your output as a BCF file (rather than as a VCF file), you can take advantage of the fact that MAPtools can also receive its input from a stream to make the BCF-to-VCF conversion:
```
bcftools view variants.bcf | ./maptools.py mbs --data D,R,Pd -m R -o mbs.txt
```

#### 5. Plotting the results
Next, run the ``plot`` command of MAPtools to automatically generate plots using the output of the previous command (i.e. the `mbs.txt` file in the example):
```
./maptools.py plot -i mbs.txt -a -o graphics_mbs/
```
For an explanation of the options available, please see the Tutorial section below.

#### 6. Optional: Re-analyzing the data with ``merge`` (recommended only for datasets with low sequencing depth and abundant polymorphisms)
If your dataset contains many polymorphic sites with low sequencing depth, re-processing the `mbs.txt` output file with the ``merge`` command can help to identify genomic regions linked to your trait of interest:
```
./maptools.py merge -i mbs.txt -o mbs_merged.txt -w 15
```
Here, `-w 15` instructs the program to combine the allele counts of 15 adjacent markers. 

Next, run the ``plot`` command again to generate plots using the file produced by ``merge`` as the input (i.e. `mbs_merged.txt`)

#### 7. Selecting a candidate region
Visually inspect the generated plots to define a broad candidate region in which your mutation is most likely to reside.

#### 8. Annotating the variants
Next, you can run the ``annotate`` command of MAPtools to automatically assess the effect of all variants present in your candidate interval. For a candidate region located on chromosome 2 between nucleotide positions 5 Mbp and 10 Mbp, you could run:
```
./maptools.py annotate -i mbs.txt -g reference_genome.gff3 -f reference_genome.fasta -R chr8:5000000-10000000 -o annotate.txt
```
Here, the `-g` flag supplies a GFF3 file containing the annotation and the `-f` flag provides a FASTA file containing the reference genome sequence. 

Please note that you cannot use the output of ``merge`` as the input for ``annotate``, as the allele information is lost in this file.

#### 9. Interpreting the results
The results of `annotate` will be saved to a text file containing detailed information on the variants identified and their putative effects. The table below corresponds to our re-analysis of a broad candidate interval on chromosome A05 of *Brassica rapa*, defined for the *nhm3-1* mutant using the data presented in Huang et al (2022). The meaning of each field (column) is explained in detail in the header of the file:
* `CHROM`: Chromosome name
* `POS`: Position in the chromosome
* `DOM`: Dominant allele
* `REC`: Recessive allele
* `DPdom_1`: Allelic depth of dominant allele in bulk 1
* `DPrec_1`: Allelic depth of recessive allele in bulk 1
* `TYPE`: Type of element in which the mutation occurs
* `ID`: ID of mutated element
* `PARENT`: The gene or mRNA containig the mutated element
* `STRAND`: Genome strand cointaining the element (+, -, .)
* `CODON_change`: Target codon before and after mutation (AGA > AAA)
* `AA_change`: Aminoacid before and after mutation (S > F)
* `INFO`: Additional information (e.g. type and effect of the variant, identity of and distance to flanking elements, etc.)

#CHROM|POS|DOM|REC|DPdom_1|DPrec_1|TYPE|ID|PARENT|STRAND|CODON_change|AA_change|INFO
---|---|---|---|---|---|---|---|---|---|---|---|---
A05|3013269|G|A|1|11|CDS|.|transcript:A05p007480.1_BraROA.1|+|GAA > AAA|E > K|effect=Nonsynonymous:missense;variant_type=substitution
A05|3160465|G|A|2|18|CDS|.|transcript:A05p007800.1_BraROA.1|+|AGA > AAA|R > K|effect=Nonsynonymous:missense;variant_type=substitution
A05|3185324|G|A|3|22|intergenic|.|.|.|.|.|right=gene:A05p007860.1_BraROA:6312;left=gene:A05p007850.1_BraROA:4983;variant_type=substitution
A05|3895591|G|A|1|13|intergenic|.|.|.|.|.|right=gene:A05p009380.1_BraROA:10963;left=gene:A05g501080.1_BraROA:4585;variant_type=substitution
A05|5135583|G|A|2|13|intron|.|transcript:A05p011720.1_BraROA.1|+|.|.|left=A05p011720.1_BraROA.1-E1:54;right=A05p011720.1_BraROA.1-E2:160;variant_type=substitution
A05|5275869|G|A|1|13|CDS|.|transcript:A05p012170.1_BraROA.1|+|GGA > GAA|G > E|effect=Nonsynonymous:missense;variant_type=substitution
A05|5311618|G|A|1|8|intergenic|.|.|.|.|.|right=gene:A05p012240.1_BraROA:2444;left=gene:A05p012230.1_BraROA:5429;variant_type=substitution
A05|5497028|C|T|2|2|intergenic|.|.|.|.|.|right=gene:A05p012640.1_BraROA:8121;left=gene:A05g501540.1_BraROA:968;variant_type=substitution
A05|5497030|C|T|2|2|intergenic|.|.|.|.|.|right=gene:A05p012640.1_BraROA:8119;left=gene:A05g501540.1_BraROA:970;variant_type=substitution
A05|5779169|G|A|0|12|intergenic|.|.|.|.|.|right=gene:A05p013290.1_BraROA:4701;left=gene:A05p013270.1_BraROA:1798;variant_type=substitution
A05|6548683|G|A|0|12|CDS|.|transcript:A05p014890.1_BraROA.1|-|ACC > ACT|T > T|effect=Synonymous:.;variant_type=substitution
A05|6687254|G|A|0|27|CDS|.|transcript:A05p015130.1_BraROA.1|+|GGA > GAA|G > E|effect=Nonsynonymous:missense;variant_type=substitution
A05|7227610|G|A|0|31|CDS|.|transcript:A05p016250.1_BraROA.1|+|GGT > GAT|G > D|effect=Nonsynonymous:missense;variant_type=substitution
A05|7360122|G|A|0|21|intergenic|.|.|.|.|.|right=gene:A05g502330.1_BraROA:17917;left=gene:A05p016520.1_BraROA:10025;variant_type=substitution
A05|7631545|G|A|3|25|intergenic|.|.|.|.|.|right=gene:A05p017070.1_BraROA:7088;left=gene:A05p017060.1_BraROA:323;variant_type=substitution
A05|7876818|G|A|0|13|intergenic|.|.|.|.|.|right=gene:A05p017450.1_BraROA:4313;left=gene:A05p017430.1_BraROA:1609;variant_type=substitution
A05|8294186|G|A|3|14|intergenic|.|.|.|.|.|right=gene:A05p018170.1_BraROA:5235;left=gene:A05p017500.1_BraROA:3564;variant_type=substitution
A05|8307063|G|A|1|10|CDS|.|transcript:A05p018190.1_BraROA.1|-|CTT > TTT|L > F|effect=Nonsynonymous:missense;variant_type=substitution
A05|9061149|G|A|1|31|intergenic|.|.|.|.|.|right=gene:A05g503000.1_BraROA:3544;left=gene:A05p019260.1_BraROA:4449;variant_type=substitution
A05|9650878|G|A|3|23|intergenic|.|.|.|.|.|right=gene:A05p020420.1_BraROA:5423;left=gene:A05p020400.1_BraROA:7195;variant_type=substitution
A05|9783762|G|A|0|24|intergenic|.|.|.|.|.|right=gene:A05p020680.1_BraROA:970;left=gene:A05p020670.1_BraROA:298;variant_type=substitution

The *nhm3-1* corresponds to the G-to-A transition mutation located at nucleotide position 5,275,869 on chromosome A05, for which all 27 reads corresponded to the recessive allele (A). See Huang *et al*. (2022) for additional details:

> Huang *et al*. (2022). *BrKAO2* mutations disrupt leafy head formation in Chinese cabbage (*Brassica rapa* L. ssp. *pekinensis*). *Theor Appl Genet* **135:** 2453-2468

##

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

In this particular experiment, each bulk has been sequenced twice. For this reason, we merge the BAM files for the ERR4463155 and ERR4463156 samples, both with white fruits:

```
samtools merge -o white.bam ERR4463155.bam ERR4463156.bam
```

```
samtools index white.bam
```

Similarly, merge the BAM files for the ERR4463153 and ERR4463154 samples, both with red fruits:

```
samtools merge -o red.bam ERR4463153.bam ERR4463154.bam
```

```
samtools index red.bam
```

#### 3. Perform the variant calling using BCFtools or GATK
The resulting VCF or BCF files of this step only need to contain the variant sites. Non-variant sites will be discarded by MAPtools later on, but larger files will increase the processing time.

3.1 With BCFtools:

Use the `mpileup` and `call` commands of BCFtools to prepare a single VCF file containing the AD (Allele Depth) field for each bulk:

```
bcftools mpileup -f Fragaria_vesca_v4.0.a1.fasta --annotate FORMAT/AD white.bam red.bam | bcftools call -mv -o variants.vcf
```

You might need to index the genome with ``samtools faidx``.


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
``-O`` specifies the output file format (default is pdf)
``-o`` specifies the output folder


#### 5. What you should get:

If you managed to follow all the steps of this tutorial, the resulting delta SNP-index plots will uncover a QTL on chromosome 1. For more information, please read the original manuscript of MAPtools, as well as the original analysis by Castillejo *et al*. (2020).


### Analysis of MBS data

For this tutorial, you will need to download the raw sequencing data from Viñegra de la Torre *et al*. (2022), which are available from the [NCBI SRA database](https://www.ncbi.nlm.nih.gov/sra) with accession numbers SRR15564670, SRR11140832, and SRR11140833. Available sequenced samples include a bulk of mutant individuals from an isogenic segregating population (SRR15564670), and two replicates of the resequencing of the non-mutagenized parent involved in the cross (SRR11140832, SRR11140833).

You will also need to download the FASTA file with the [reference genome sequence of *Arabis alpina*](http://www.arabis-alpina.org/data/ArabisAlpina/assemblies/V5.1/), and the [genome annotation file in GFF3 format](http://www.arabis-alpina.org/data/ArabisAlpina/assemblies/V5.1/)

#### 1. Map the reads to the reference genome using your read mapper of choice

1.1. With Bowtie2:

Index the reference genome using bowtie2-build:

```
bowtie2-build Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz INDEX
```

Map the paired-end reads to the reference genome using bowtie2. Repeat this step for each available sample:

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

#### 4. What you should get:
In this tutorial, we have analyzed one (out of 3) allelic mutants described in the original paper. If you followed all the steps of this tutorial, the plots should uncover a bias in the allele frequencies of chromosome 8. If you analyze the reads for the other 2 mutants described in the original paper, you will find lesions in the same gene! For more information, please read our MAPtools paper and the original study by Viñegra de la Torre *et al*. (2022).


## **Citation**

If you find this software useful, please cite:

César Martínez-Guardiola, Ricardo Parreño, Héctor Candela (2024).
MAPtools: command-line tools for mapping-by-sequencing and QTL-Seq analysis and visualization. *Plant Methods* **20,** 107 (2024). 
https://doi.org/10.1186/s13007-024-01222-2
