
fields = {'#CHROM': str, 'POS': int, 'DOM': str, 'REC': str, 'DPdom_1': int, 'DPrec_1': int, 'DPdom_2': int, 'DPrec_2': int, 'SNPidx1': float,
          'SNPidx2': float, 'DELTA': float, 'MAX_SNPidx2': float, 'PVALUE': float, 'log10PVALUE': float, 'FISHER': float, 'BOOST': float, 'ED': float, 'G':float, 'CI95': float}
RANG = 100
#MUT = 'EMS'
blind_mbs = {'dots':('#000000','black'), 'MAX_SNPidx2':('#f781bf','pink'),'log10PVALUE':('#f0e442','yellow'),'SNPidx1':('#00ffff','light blue'),'SNPidx2':('#ff9a31','orange'),'BOOST':('#00ffff','light blue')}
norm_mbs = {'dots':('royalblue','blue'), 'MAX_SNPidx2':('#f781bf','pink'),'log10PVALUE':('tomato','red'),'SNPidx1':('limegreen','green'),'SNPidx2':('#FF33F9' ,'pink'),'BOOST':('limegreen','green')}

palettes_mbs = {'1':blind_mbs,'2':norm_mbs}
#TODO - Configure qtl palettes
norm_qtl = {'dots':('royalblue','blue'), 'MAX_SNPidx2':('limegreen','green')
,'log10PVALUE':('tomato','red'),'SNPidx1':('limegreen','green'),'SNPidx2':('#FAC205','yellow'),
'DELTA':('#f781bf','pink'), 'mvg':('tomato','red'), 'ci':('#ff7f00','orange')
}

blind_qtl={'dots':('#000000','black'), 'MAX_SNPidx2':('#f781bf','pink')
,'log10PVALUE':('#f0e442','yellow'),'SNPidx1':('#00ffff','light blue'),'SNPidx2':('#ff9a31','orange'),
'DELTA':('#f781bf','pink'), 'mvg':('#f0e442','yellow'), 'ci':('#f781bf','pink')
}
palettes_qtl = {'1':blind_qtl,'2':norm_qtl}

titles_mbs = {1:('p-values of two-tailed Fisher\'s exact tests along the chromosomes ({}).','pvalues_multiH' ),
          2:('p-values of two-tailed Fisher\'s exact tests in chromosome {}.', 'pvalue_{}'),
          3:('Frequency of alleles in the pool of plants exhibiting the wild-type phenotype, in chromosome {}.', 'AF1&2_{}', 'AF1_{}'),
          4:('Frequency of alleles in the pool of plants exhibiting the mutant phenotype, in chromosome {}.', 'AF2&1_{}', 'AF2_{}'),
          5:('Maximum of allele frequency of plants exhibiting the mutant phenotype, in chromosome {}.', 'MAXAF2_{}'),
          6:('p-values of two-tailed Fisher\'s exact tests.','pvalues_multiV'),
          7:('Frequency of alleles in the pool of plants exhibiting the wild-type phenotype.', 'AF1&2_multiV','AF1_multiV'),
          8:('Frequency of alleles in the pool of plants exhibiting the mutant phenotype.', 'AF2&1_multiV','AF2_multiV'),
          9:('Maximum frequency of alleles in the pool of plants exhibiting the mutant phenotype.', 'MAXAF2_multiV'),
          10:('Mapping-by-sequencing of chromosome {}.', 'AF1&2&pval_{}'),
          11:('Single nucleotide polymorphism (SNP)-index plots of chromosome {}','SNPidx_{}')
}

titles_qtl = {1:('p-values of two-tailed Fisher\'s exact tests along the chromosomes ({}).','pvalues_multiH' ),
          2:('p-values of two-tailed Fisher\'s exact tests in chromosome {}.', 'pvalue_{}'),
          3:('Frequency of alleles in the pool of plants exhibiting the high phenotype, in chromosome {}.', 'AF1&2_{}', 'AF1_{}'),
          4:('Frequency of alleles in the pool of plants exhibiting the low phenotype, in chromosome {}.', 'AF2&1_{}', 'AF2_{}'),
          5:('$\Delta$ SNP-index for chromosome {}.','Delta_{}'),
          6:('p-values of two-tailed Fisher\'s exact tests.','pvalues_multiV'),
          7:('Frequency of alleles in the pool of plants exhibiting the high phenotype.', 'AF1&2_multiV','AF1_multiV'),
          8:('Frequency of alleles in the pool of plants exhibiting the low phenotype.', 'AF2&1_multiV','AF2_multiV'),
          9:('$\Delta$ SNP-index plot with data from both pools.','Delta_multiV'),
          11:('Single nucleotide polymorphism (SNP)-index plots of chromosome {}','SNPidx_{}'),
          10:('Euclidean Distance between frequencies in both pools', 'ED_{}','ED_multiV'),
          12:('G-statistic of the two-phase sampling inherent to both pools, calculated using the procedure used in (Magwene, et al. 2011) .','G_{}','G_multiV'),
          13:('QTL-Seq comparison plot of chromosome {}.', 'QTL-Seq_{}')
}

lines_mbs = {-3:'Each {} dot indicates the  maximum allele frequency of a biallelic SNP segregating in the population, as determined for the pool of mutants.',
        -2:'Each {} dot indicates the allele frequency of a biallelic SNP segregating in the population, as determined for the pool of wild-type sibilings.',
        -1:'Each {} dot indicates the allele frequency of a biallelic SNP segregating in the population, as determined for the pool of mutants.',
         0:'Each {} dot indicates the p-value of a biallelic SNP segregating in the population, as determined using the data from both pools.',
         1:'The {} line indicates the weighted moving average of the p-values of {} adjacent SNPs.',
	     2:'The dashed line marks the significance threshold calculated using the Bonferroni correction, considering n={} chromosomal locations have been tested.',
         3:'The {} line indicates the moving average of the allele frequencies at {} adjacent SNPs, in a pool of plants exhibiting the wild-type phenotype.',
         4:'The {} line indicates the moving average of the allele frequencies at {} adjacent SNPs, in a pool of plants exhibiting the mutant phenotype.',
         5:'The {} line indicates the moving average of the maximum allele frequencies at {} adjacent SNPs, in a pool of plants exhibiting the mutant phenotype.',
         6:'The {} dashed line is the moving average of Boost at {} adjacent SNPs. Boost was calculated from the allele frecuency values of the recessive pool, using the procedure described in SHOREmap (Schneeberger, et al. 2009).',
         7:'(a) Frequency of haplotypes inherited from the mutant parent exhibiting the wild-type phenotype.'\
            '(b) Frequency of haplotypes inherited from the mutant parent exhibiting the mutant phenotype.'\
                '(c) p-values of two-tailed Fisher\'s exact tests.'\
                    'Each {} dot indicates the allele frequency (a,b) or the p-value (c) of biallelic SNP segregating population.'
}

lines_qtl = {-3:'Each {} dot was calculated substracting the SNP-index value of the low pool from the high pool.',
        -2:'Each {} dot indicates the allele frequency of a biallelic SNP segregating in the population, as determined for the pool of high phenotype.',
        -1:'Each {} dot indicates the allele frequency of a biallelic SNP segregating in the population, as determined for the pool of low phenotype.',
         0:'Each {} dot indicates the p-value of a biallelic SNP segregating in the population, as determined using the data from both pools.',
         1:'The {} line indicates the weighted moving average of the p-values of {} adjacent SNPs.',
	      2:'The dashed line marks the significance threshold calculated using the Bonferroni correction, considering n={} chromosomal locations have been tested.',
         3:'The {} line indicates the moving average of the allele frequencies at {} adjacent SNPs, in a pool of plants exhibiting the high phenotype.',
         4:'The {} line indicates the moving average of the allele frequencies at {} adjacent SNPs, in a pool of plants exhibiting the low phenotype.',
         5:'The {} line indicates the moving average of the $\Delta$ SNP-index at {} adjacent SNPs.',
         6:'(a) SNP-index plot of the pool of high phenotype plants. (b) SNP-index plot of the pool of low phenotype plants'\
            ' (c) $\Delta$ SNP-index plot, substracting the SNP-index value of the low pool from the high pool.',
         7:'The {} region shows the confidence interval for a proportion difference at a significance level of 5%, adjusted by the Bonferroni correction (n={}).',
         8:'The {} line indicates the moving average for each of the plotted variables, at {} adjacent SNPs.',
         9:'Each {} dot represents the Euclidean distance for each individual marker',
         10:'The {} line is the smoothed statistic ED100^4 for Euclidean distance, calculated by summing of 100 individual markers raised to the fourth power as described in (de la Fuente Cantó, et al. 2022).',
         11:'Each {} dot represents the G-statistic for each individual marker',
         12:'The {} line indicates the weighted moving average of the G-statistic at {} adjacent SNPs.',
         13:'(a) G-statistic for each individual SNP. (b) Euclidean distance for each individual SNP and the smoothed statistic ED100^4.'\
            '(c) $\Delta$ SNP-index plot, substracting the SNP-index value of the low pool from the high pool.'\
               '(d) p-values of two-tailed Fisher\'s exact tests.'
        
         
}
v_maptools='Maptools version: 1.00'
v_annotate= 'Annotate variants version: 1.00'
v_mbs= 'Mapping by Sequencing analysis version: 1.00'
v_mbsplot = 'Mapping by Sequencing plots version: 1.00'
v_qtl= 'QTL-Seq analysis version: 1.00'
v_qtlplot= 'QTL-Seq plots version: 1.00'
v_merge='Merge mbs or qtl analysis version: 1.00'

variable_descriptions={
   '#CHROM':'##CHROM=<ID=*, Description=\"Chromosome name.\">\n',\
   'POS':'##POS=<ID=*, Number=1, Type=Integer, Description=\"Position in the chromosome.\">\n',\
   'DOM':'##DOM=<ID=*, Description=\"Dominant allele.\">\n',\
   'REF':'##REF=<ID=*, Description=\"Reference allele.\">\n',\
   'REC':'##REC=<ID=*, Description=\"Recessive allele.\">\n',\
   'ALT':'##ALT=<ID=*, Description=\"Represents allele, other than observed.\">\n',\
   'DPdom_1':'##DPdom_1=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of dominant allele in pool 1.\">\n',\
   'DPrec_1':'##DPrec_1=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of recessive allele in pool 1.\">\n',\
   'DPdom_2':'##DPdom_2=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of dominant allele in pool 2.\">\n',\
   'DPrec_2':'##DPrec_2=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of recessive allele in pool 2.\">\n',\
   'DPref_1':'##DPref_1=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of reference allele in pool 1.\">\n',\
   'DPalt_1':'##DPalt_1=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of alternative allele in pool 1.\">\n',\
   'DPref_2':'##DPref_2=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of reference allele in pool 2.\">\n',\
   'DPalt_2':'##DPalt_2=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of alternative allele in pool 2.\">\n',\
   'SNPidx1':'##SNPidx1=<ID=*, Number=1, Type=Float, Description=\"Allele frequency of pool 1.\">\n',\
   'SNPidx2':'##SNPidx2=<ID=*, Number=1, Type=Float, Description=\"Allele frequency of pool 1.\">\n',\
   'MAX_SNPidx2':'##MAX_SNPidx2=<ID=*, Number=1, Type=Float, Description=\"Maximum allele frecuency ratio of target pool.\">\n',\
   'FISHER':'##FISHER=<ID=*, Number=1, Type=Float, Description=\"Two-tailed Fisher exact test.\">\n',\
   'BOOST':'##BOOST=<ID=*, Number=1, Type=Float, Description=\"Scarcity-values, calculated from the allele frecuency values of the recessive pool.\">\n',\
   'PVALUE':'##PVALUE=<ID=*, Number=1, Type=Float, Description=\"P-values of Fisher test.\">\n',\
   'log10PVALUE':'##log10PVALUE=<ID=*, Number=1, Type=Float, Description=\"Base 10 logarithm of p-value.\">\n',\
   'DELTA':'##DELTA=<ID=*, Number=1, Type=Float, Description=\"Difference between SNPidx1 minus SNPidx2 (1 to -1).\">\n',\
   'ED':'##ED=<ID=*, Number=1, Type=Float, Description=\"Euclidean distance.\">\n',\
   'G':'##G=<ID=*, Number=1, Type=Float, Description=\"G-statistic.\">\n',\
   'TYPE':'##TYPE=<ID=*, Number=0, Type=Flag, Description=\"Type of element in which the mutation occurs.\">\n',\
   'ID':'##ID=<ID=*, Number=0, Type=Flag, Description=\"ID of mutated element.\">\n',\
   'PARENT':'##PARENT=<ID=*, Number=0, Type=Flag, Description=\"The gene or mRNA containig the mutated element.\">\n',\
   'STRAND':'##STRAND=<ID=*, Number=0, Type=Flag, Description=\"Genome strand cointaining the element (+, -, .).\">\n',\
   'PHASE':'##PHASE=<ID=*, Number=0, Type=Integer, Description=\"CDS phase (0,1,2)\".>\n',\
   'CODON_change':'##CODON_change=<ID=*, Number=0, Type=Flag, Description=\"Target codon before and after mutation (AGA > AAA).\">\n',\
   'AA_change':'##AA_change=<ID=*, Number=0, Type=Flag, Description=\"Aminoacid before and after mutation (S > F).\">\n',\
   'INFO':'##INFO=<ID=effect, Number=0, Type=Flag, Description=\"Type of aminoacid change effect in exons (Nonsynonymous or Synonymous:missense or nonsense or nonstop), in 5\'UTR (new_ATG:in_frame or out_of_frame:truncated_protein or elongated_protein) or if affects to non-coding gene (non_coding:exonic or intronic).\">\n'\
          '##INFO=<ID=left, Number=0, Type=Flag, Description=\"Gene or Exon immediately left (towards lower POS) and distance on bp(element:distance).\">\n'\
          '##INFO=<ID=right, Number=0, Type=Flag, Description=\"Gene or Exon immediately right (towards higher POS) and distance on bp(element:distance).\">\n'\
          '##INFO=<ID=5_splice_site, Number=0, Type=Flag, Description=\"The mutation is 3 or less bp from the upstream splicing site (exon or intron _boundary:ID_CDS).">\n'\
          '##INFO=<ID=3_splice_site, Number=0, Type=Flag, Description=\"The mutation is 3 or less bp from the downstream splicing site (exon or intron _boundary:ID_CDS).">\n'\
   
}
