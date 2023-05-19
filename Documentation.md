# Maptools Documentation Page
## VERSION
This documentation page was last updated **2023-05-19** an refers to MAPtools version **1.00**.

## LIST OF COMMANDS
MAPtools consists of several commands that are used in a similar way to the commands in samtools and bcftools.

To see the list of available commands, run maptools with no arguments. To see the options for each command run maptools *command* without options. The current version includes the following commands:

* ``mbs`` : analyzes results of mapping-by-sequencing.
* ``mbsplot`` : automatically generates various types of figures (and their legends) from the results of ``mbs``.
* ``qtl`` : analyzes results of QTL-seq.
* ``qtlplot`` : generates figures from the results of ``qtl``, similar to ``mbsplot``.
* ``merge``: reanalyzes data using pseudomarkers created by combining the results of multiple individual markers (similar to counting haplotype alleles instead of marker alleles). The results of ``merge`` can also be represented with ``mbsplot`` and ``qtlplot``, depending on whether they are of type ``mbs`` or ``qtl``.
* ``annotate``: evaluates the effect of SNPs and indels contained within a specified interval.

## COMMANDS AND OPTIONS
### maptools.py mbs [options]
**-i, --input** FILE
 - Requiere un archivo vcf como entrada, procedente de un variant calling. Esta opción puede obviarse si el input viene de un pipe.

**

    
