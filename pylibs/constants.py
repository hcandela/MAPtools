blind_mbs = {'dots':('#000000','black'), 'MAX_SNPidx2':('#f781bf','pink'),'log10PVALUE':(),'SNPidx1':('#00ffff','light blue'),'SNPidx2':('#ff9a31','orange'),'BOOST':('#00ffff','light blue')}
norm_mbs = {'dots':('royalblue','blue'), 'MAX_SNPidx2':('#f781bf','pink'),'log10PVALUE':('tomato','red'),'SNPidx1':('limegreen','green'),'SNPidx2':('#FF33F9' ,'pink'),'BOOST':('limegreen','green')}

palettes_mbs = {'1':blind_mbs,'2':norm_mbs}

norm_qtl = {'dots':('royalblue','blue'), 'MAX_SNPidx2':('limegreen','green')
,'log10PVALUE':('tomato','red'),'SNPidx1':('limegreen','green'),'SNPidx2':('#FAC205','yellow'),
'DELTA':('#f781bf','pink'), 'mvg':('tomato','red'), 'ci':('#ff7f00','orange')
}

palettes_qtl = {'2':norm_qtl}

titles_mbs = {1:('p-values of two-tailed Fisher\'s exact tests along the chromosomes ({}).','pvalues_multiH' ),
          2:('p-values of two-tailed Fisher\'s exact tests in chromosome {}.', 'pvalue_chr{}'),
          3:('Frequency of alleles in the pool of plants exhibiting the wild-type phenotype, in chromosome {}.', 'AF1&2_chr{}', 'AF1_chr{}'),
          4:('Frequency of alleles in the pool of plants exhibiting the mutant phenotype, in chromosome {}.', 'AF2&1_chr{}', 'AF2_chr{}'),
          5:('Maximum of allele frequency of plants exhibiting the mutant phenotype, in chromosome {}.', 'MAXAF2_chr{}'),
          6:('p-values of two-tailed Fisher\'s exact tests.','pvalues_multiV'),
          7:('Frequency of alleles in the pool of plants exhibiting the wild-type phenotype.', 'AF1&2_multiV','AF1_multiV'),
          8:('Frequency of alleles in the pool of plants exhibiting the mutant phenotype.', 'AF2&1_multiV','AF2_multiV'),
          9:('Maximum frequency of alleles in the pool of plants exhibiting the mutant phenotype.', 'MAXAF2_multiV'),
          10:('Mapping-by-sequencing of chromosome {}.', 'AF1&2&pval_chr{}'),
          11:('Single nucleotide polymorphism (SNP)-index plots of chromosome {}','SNPidx_chr{}')
}

titles_qtl = {1:('p-values of two-tailed Fisher\'s exact tests along the chromosomes ({}).','pvalues_multiH' ),
          2:('p-values of two-tailed Fisher\'s exact tests in chromosome {}.', 'pvalue_chr{}'),
          3:('Frequency of alleles in the pool of plants exhibiting the high phenotype, in chromosome {}.', 'AF1&2_chr{}', 'AF1_chr{}'),
          4:('Frequency of alleles in the pool of plants exhibiting the low phenotype, in chromosome {}.', 'AF2&1_chr{}', 'AF2_chr{}'),
          5:('$\Delta$ SNP-index for chromosome {}.','Delta_chr{}'),
          6:('p-values of two-tailed Fisher\'s exact tests.','pvalues_multiV'),
          7:('Frequency of alleles in the pool of plants exhibiting the high phenotype.', 'AF1&2_multiV','AF1_multiV'),
          8:('Frequency of alleles in the pool of plants exhibiting the low phenotype.', 'AF2&1_multiV','AF2_multiV'),
          9:('$\Delta$ SNP-index plot with data from both pools.','Delta_multiV'),
          11:('Single nucleotide polymorphism (SNP)-index plots of chromosome {}','SNPidx_chr{}')
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
         7:'(a) Frequency of haplotypes inherited from the mutant parent exhibiting the wild-type phenotype.\n'\
            '(b) Frequency of haplotypes inherited from the mutant parent exhibiting the mutant phenotype.\n'\
                '(c) p-values of two-tailed Fisher\'s exact tests.\n'\
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
         8:'The {} line indicates the moving average for each of the plotted variables, at {} adjacent SNPs.'
        
         
}
