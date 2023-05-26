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


In this experiment, each bulk has been sequenced twice. For this reason, we need to merge the BAM files for the ERR4463155 and ERR4463156 samples, both of which with white fruits:

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

* First, we analyze the data with ``qtl`` option, assuming that the white and red are the extreme phenotypes "low" and "high", respectively. ``-I`` option is used to ignore indels for this example:

```
./maptools qtl --data L,H -i variant_calling_out.vcf -o qtl_output.txt -I
```

* Now we can plot the results using ``qtlplot``:

```
./maptools qtlplot -i qtl_output.txt --captions -m -c Fvb1,Fvb2,Fvb3,Fvb4,Fvb5,Fvb6,Fvb7 -A 800 -a --bonferroni --ci95 -o ../graphics_qtl/
```

* Explanation of options:

``--captions`` Generates figure captions.
``-m`` Multi-chromosome plots.
``-c`` List of selected chromosomes for visualization.
``-A`` Number of adjacent markers to calculate moving average lines.
``-a`` Generates all possible plots.
``--bonferroni`` Show bonferroni test line in p-value plots.
``--ci95`` Show the 95% confidence interval in delta plots.
``-o`` Output folder.

#### Analysis of MBS data

For this tutorial, you will need to download the raw sequencing data from Viñegra de la Torre *et al*. (2022), which are available from the [NCBI SRA database](https://www.ncbi.nlm.nih.gov/sr) with accession numbers SRR15564670, SRR11140832, and SRR11140833. You will also need to download the FASTA file with the [reference genome sequence of *Arabis alpina*](http://www.arabis-alpina.org/data/ArabisAlpina/assemblies/V5.1/), and the [genome annotation file in GFF3 format](http://www.arabis-alpina.org/data/ArabisAlpina/assemblies/V5.1/)

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
```


In this experiment, the parental sample has been sequenced twice. For this reason, we need to merge the BAM files for the SRR11140832 and SRR11140833 files:

```
samtools merge -o parental.bam SRR11140832.bam SRR11140833.bam
```

```
samtools index parental.bam
```

```
samtools index mutant.bam
```

#### 3. Perform the variant calling using BCFtools or GATK

3.1 With BCFtools:

Use the `mpileup` and `call` commands of BCFtools to prepare a single VCF file containing the AD (Allele Depth) field for each bulk:

```
bcftools mpileup -f Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz  --annotate FORMAT/AD mutant.bam parental.bam | bcftools call -mv -o variants.vcf
```

You might need to index the genome with ``samtools faidx``

3.2 With GATK:

If you plan to do the variant calling with GATK, note that GATK requires BAM files with the RG field set. You can add this field with the `addreplacerg` command of samtools:

```
samtools addreplacerg -r ID:mutant -r SM:mutant -o fixed_mutant.bam mutant.bam

samtools addreplacerg -r ID:parental -r SM:parental -o fixed_parental.bam parental.bam
```

Use the `HaplotypeCaller` command of GATK to prepare a single VCF file containing the AD (Allele Depth) field for each bulk:

```
gatk HaplotypeCaller --output variants.vcf -I fixed_mutant.bam -I fixed_parental.bam --reference Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz

```

#### 4. Running MAPtools:

* First, we analyze the data with ``mbs`` option assuming that the mutant phenotype is in the recesive pool. Also, we have processed the samples of the recesive pool and the parental dominant resequencing with BCFtools, in that specific order:

```
./maptools.py mbs --data R,Pd -m R --EMS --parental-filter -o {mbs_out.txt} -i {variant_calling_mbs.vcf}
```

* Explanation of options:

``--data``: defines the type of data. In this case the ``bcftools mpileup`` has processed the samples in this order: recesive sample pool (R), and parental dominant (Pd).
``-m R``: indicates which pool carries the mutant phenotype.
``--EMS``: filter out SNPs not caused by EMS.
``--parental-filter``: filter out variants present in the parental sequencing.
``-o``: output file.

#### 6. visualizing data with MAPtools:

We generate a plot of the data generated with ``mbsplot``. In this case it is recommended to specify the chromosomes that are going to be ploted, as there are a lot of individual contigs in the reference fasta.

```
./maptools.py mbsplot -i {mbs_out.txt} -m --captions -c chr1,chr2,chr3,chr4,chr5,chr6,chr6,chr7,chr8 -A 10 -b 20 -a -o {../graphics_mbs/}
```

* Explanation of options:

``-m``: multi-chromosome plots.
``--captions``: generates figure captions.
``-c``: list of selected chromosomes for visualization.
``-A``: number of adjacent markers to calculate moving average lines.
``-b``: generates a boost plot with an moving average line of the given INT (in this case, 20).
``-a``: generates all possible plots.
``-o``: output directory.

#### 7. Analyze MBS data with MAPtools:

As we can observe in the plots, the mutation is located in chromosome 8 and approximately between positions 14250000 and 16700000, so we can analyze this region:

```
./maptools.py annotate --data R,Pd -m R --EMS --parental-filter -i {variant_calling_mbs.vcf} -g {gff_arabis.gff3} -f Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz -R chr8:14250000-16700000 -o {annotate_mbs.txt}
```

* Explanation of options:

``--data R,Pd``: defines the type of data. In this case, the samples are the the recessive pool (R) and the parental dominant (Pd).
``-m R``: indicates which pool carries the mutant phenotype. It is not necessary in this example, but it is present to improve clarity.
``--EMS``: filter out SNPs not caused by EMS.
``--parental-filter``: filter out variants present in the parental sequencing.
``-i``: input file in VCF format'
``-g``: GFF annotation file.
``-f``: reference genome file.
``-R``: region where the analysis will be performed.
``-o``: output file.

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
