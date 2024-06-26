v_maptools='MAPtools version 1.00'
v_annotate= 'MAPtools annotate command version 1.00'
v_mbs= 'MAPtools mbs command version 1.00'
v_plot = 'MAPtools plot command version 1.00'
v_qtl= 'MAPtools qtl command version 1.00'
v_merge='MAPtools merge command version 1.00'

fields = {'#CHROM': str, 'POS': int, 'DOM': str, 'REC': str, 'DPdom_1': int, 'DPrec_1': int, 'DPdom_2': int, 'DPrec_2': int, 'SNPidx1': float,
          'SNPidx2': float, 'DELTA': float, 'MAX_SNPidx2': float, 'PVALUE': float, 'log10PVALUE': float, 'FISHER': float, 'BOOST': float, 'ED': float, 'G':float, 'CI95': float}
RANG = 100

TITLES = {#AF MBS
          0:('Frequency of alleles in the bulk of individuals exhibiting the dominant phenotype, in chromosome {}.', 'AF1&2_{}', 'AF1_{}'),
          1:('Frequency of alleles in the bulk of individuals exhibiting the recessive phenotype, in chromosome {}.', 'AF2&1_{}', 'AF2_{}'),
          2:('Frequency of alleles in the bulk of individuals exhibiting the dominant phenotype.', 'AF1&2_multiV','AF1_multiV', 'AF1_Manhattan'),
          3:('Frequency of alleles in the bulk of individuals exhibiting the recessive phenotype.', 'AF2&1_multiV','AF2_multiV','AF2_Manhattan'),
          4:('Frequency of alleles in both bulks.', 'AF2&1_Manhattan'),
          #AF QTL
          5:('Frequency of alleles in the bulk of individuals exhibiting the high phenotype, in chromosome {}.', 'AF1&2_{}', 'AF1_{}'),
          6:('Frequency of alleles in the bulk of individuals exhibiting the low phenotype, in chromosome {}.', 'AF2&1_{}', 'AF2_{}'),
          7:('Frequency of alleles in the bulk of individuals exhibiting the high phenotype.', 'AF1&2_multiV','AF1_multiV', 'AF1_Manhattan'),
          8:('Frequency of alleles in the bulk of individuals exhibiting the low phenotype.', 'AF2&1_multiV','AF2_multiV','AF2_Manhattan'),
          #DELTA
          9:('\u0394 SNP-index for chromosome {}.','delta_{}'),
          10:('\u0394 SNP-index plot with data from both bulks.','delta_multiV','delta_Manhattan'),
          #p-value
          11:('p-values of two-tailed Fisher\'s exact tests in chromosome {}.', 'pvalue_{}'),
          12:('p-values of two-tailed Fisher\'s exact tests.', 'pvalue_multiV','pvalue_Manhattan'),
          #ED
          13:('Euclidean Distance between frequencies in both bulks in chromosome {}.', 'ED_{}'),
          14:('Euclidean Distance between frequencies in both bulks.','ED_multiV', 'ED_Manhattan'),
          #G-statistic
          15:('G-statistic of the two-phase sampling inherent to both bulks in chromosome {}, calculated using the procedure used in Magwene et al. (2011).','G_{}'),
          16:('G-statistic of the two-phase sampling inherent to both bulks, calculated using the procedure used in Magwene et al. (2011).','G_multiV', 'G_Manhattan'),
          #Combined
          17:('Different statistics place on chromosome {}.', 'Combined_statistics_{}'),
          #Maximum
          18:('Maximum of allele frequency of individuals exhibiting the recessive phenotype, in chromosome {}.', 'MAXAF2_{}'),
          19:('Maximum of allele frequency of individuals exhibiting the recessive phenotype', 'MAXAF2_multiV', 'MAXAF2_Manhattan')
}
         
LINES = {#AF MBS
         0:'Each {} dot indicates the allele frequency of a biallelic polymorphism segregating in the population, as determined for the bulk exhibiting the dominant phenotype.',
         1:'Each {} dot indicates the allele frequency of a biallelic polymorphism segregating in the population, as determined for the bulk exhibiting the recessive phenotype.',
         2:'The {} line represents the weighted average of the allele frequencies of markers within a {}-bp range, in a bulk of individuals exhibiting the dominant phenotype.',
         3:'The {} line indicates the moving average of the allele frequencies at {} adjacent polymorphisms, in a bulk of individuals exhibiting the dominant phenotype.',
         4:'The {} line represents the weighted average of the allele frequencies of markers within a {}-bp range, in a bulk of individuals exhibiting the recessive phenotype.',
         5:'The {} line indicates the moving average of the allele frequencies at {} adjacent polymorphisms, in a bulk of individuals exhibiting the recessive phenotype.',
         #AF QTL
         6:'Each {} dot indicates the allele frequency of a biallelic polymorphism segregating in the population, as determined for the bulk exhibiting the high phenotype.',
         7:'Each {} dot indicates the allele frequency of a biallelic polymorphism segregating in the population, as determined for the bulk exhibiting the low phenotype.',
         8:'The {} line represents the weighted average of the allele frequencies of markers within a {}-bp range, in a bulk of individuals exhibiting the high phenotype.',
         9:'The {} line indicates the moving average of the allele frequencies at {} adjacent polymorphisms, in a bulk of individuals exhibiting the high phenotype.',
         10:'The {} line represents the weighted average of the allele frequencies of markers within a {}-bp range, in a bulk of individuals exhibiting the low phenotype.',
         11:'The {} line indicates the moving average of the allele frequencies at {} adjacent polymorphisms, in a bulk of individuals exhibiting the low phenotype.',
         #p-value
         12:'Each {} dot indicates the {}log(p-value) of a biallelic marker segregating in the population, as determined using the data from both bulks.',
         13:'The {} line represents the weighted average of the {}log(p-values) of markers within a {}-bp range.',
         14:'The {} line indicates the weighted moving average of the {}log(p-values) of {} adjacent markers.',
         15:'The dashed line marks the significance threshold calculated using the Bonferroni correction, considering that n={} chromosomal locations have been tested.',
         #delta
            #MBS
         16:'Each {} dot was calculated subtracting the allele frequency value of the recessive bulk from the dominant bulk.',
         17:'Each {} dot was calculated subtracting the allele frequency value of the recessive bulk from the dominant bulk, and obtaining its absolute value.',
            #QTL
         18:'Each {} dot was calculated subtracting the allele frequency value of the low bulk from the high bulk.',
         19:'Each {} dot was calculated subtracting the allele frequency value of the low bulk from the high bulk, and obtaining its absolute value.',
         20:'The {} line indicates the moving average of the \u0394 SNP-index at {} adjacent polymorphisms.',
         21:'The {} line represents the weighted average of the \u0394 SNP-index of markers within a {}-bp range.',
         22:'The shaded region shows the confidence interval for a proportion difference at a significance level of 5%, adjusted by the Bonferroni correction (n={}).',
         #ED
         23:'Each {} dot represents the Euclidean distance for each individual marker.',
         24:'The {} line is the smoothed statistic ED100^4 for Euclidean distance, calculated by summing of 100 individual markers raised to the fourth power as described in de la Fuente Cantó, et al. (2022).',
         #G
         25:'Each {} dot represents the G-statistic for each individual marker.',
         26:'The {} line represents the weighted average of the G-statistic of markers within a {}-bp range.',
         27:'The {} line indicates the weighted moving average of the G-statistic at {} adjacent polymorphisms.',
         #Combined
         28:' Unless otherwise stated, continuous lines represents the weighted moving average calculated using the markers within a {}-bp range.',
         29:' Unless otherwise stated, continuous lines represents the weighted moving average calculated using a sliding window of {} adjacent markers.',
            #MBS
         30:'Each dot corresponds to an individual biallelic marker segregating in the population.'\
            '{}'\
            '(a) SNP-index plot for the bulk of dominant phenotype individuals. (b) SNP-index plot for the bulk of recessive phenotype individuals.'\
            ' (c) \u0394 SNP-index plot, subtracting the SNP-index value of the recessive bulk from that of the dominant bulk{}. {}'\
            ' (d) Euclidean distance for each individual marker (dots). {}'\
            ' (e) G-statistic for each individual marker.'\
            ' (f) p-values of two-tailed Fisher\'s exact tests. {}',
            #QTL
         31:'Each dot corresponds to an individual biallelic marker segregating in the population.'\
            '{}'\
            '(a) SNP-index plot for the bulk of high phenotype individuals. (b) SNP-index plot for the bulk of low phenotype individuals.'\
            ' (c) \u0394 SNP-index plot, subtracting the SNP-index value of the low bulk from that of the high bulk{}. {}'\
            ' (d) Euclidean distance for each individual marker (dots). {}'\
            ' (e) G-statistic for each individual marker.'\
            ' (f) p-values of two-tailed Fisher\'s exact tests. {}',         
         #MAX
         32: 'Each {} dot indicates the  maximum allele frequency of a biallelic polymorphism segregating in the population, as determined for the bulk with the recessive phenotype.',
         33: 'The {} line indicates the moving average of the maximum allele frequencies at {} adjacent polymorphisms, in a bulk of individuals exhibiting the recessive phenotype.',
         34: 'The {} line represents the weighted average of the maximum allele frequencies of markers within a {}-bp range, in a bulk of individuals exhibiting the recessive phenotype.',
         #BOOST
         35:'The {} dashed line is the weighted average of Boost of markers within a  {}-bp range. Boost was calculated from the allele frecuency values of the recessive bulk, using the procedure described in SHOREmap (Schneeberger, et al. 2009).',
         36:'The {} dashed line is the moving average of Boost at {} adjacent polymorphisms. Boost was calculated from the allele frecuency values of the recessive bulk, using the procedure described in SHOREmap (Schneeberger, et al. 2009).'

}

variable_descriptions={
   '#CHROM':'##CHROM=<ID=*, Description=\"Chromosome name.\">\n',\
   'POS':'##POS=<ID=*, Number=1, Type=Integer, Description=\"Position in the chromosome.\">\n',\
   'DOM':'##DOM=<ID=*, Description=\"Dominant allele.\">\n',\
   'REF':'##REF=<ID=*, Description=\"Reference allele.\">\n',\
   'REC':'##REC=<ID=*, Description=\"Recessive allele.\">\n',\
   'ALT':'##ALT=<ID=*, Description=\"Represents allele, other than observed.\">\n',\
   'HIGH':'##HIGH=<ID=*, Description=\"Allele from the high bulk.\">\n',\
   'LOW':'##LOW=<ID=*, Description=\"Allele from the low bulk.\">\n',\
   'DPdom_1':'##DPdom_1=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of dominant allele in bulk 1.\">\n',\
   'DPrec_1':'##DPrec_1=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of recessive allele in bulk 1.\">\n',\
   'DPdom_2':'##DPdom_2=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of dominant allele in bulk 2.\">\n',\
   'DPrec_2':'##DPrec_2=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of recessive allele in bulk 2.\">\n',\
   'DPref_1':'##DPref_1=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of reference allele in bulk 1.\">\n',\
   'DPalt_1':'##DPalt_1=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of alternative allele in bulk 1.\">\n',\
   'DPref_2':'##DPref_2=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of reference allele in bulk 2.\">\n',\
   'DPalt_2':'##DPalt_2=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of alternative allele in bulk 2.\">\n',\
   'DPhigh_1':'##DPhigh_1=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of high allele in bulk 1.\">\n',\
   'DPlow_1':'##DPlow_1=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of low allele in bulk 1.\">\n',\
   'DPhigh_2':'##DPhigh_2=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of high allele in bulk 2.\">\n',\
   'DPlow_2':'##DPalt_2=<ID=*, Number=1, Type=Integer, Description=\"Allelic depth of low allele in bulk 2.\">\n',\
   'SNPidx1':'##SNPidx1=<ID=*, Number=1, Type=Float, Description=\"Allele frequency of bulk 1.\">\n',\
   'SNPidx2':'##SNPidx2=<ID=*, Number=1, Type=Float, Description=\"Allele frequency of bulk 2.\">\n',\
   'MAX_SNPidx2':'##MAX_SNPidx2=<ID=*, Number=1, Type=Float, Description=\"Maximum allele frecuency ratio of target bulk.\">\n',\
   'FISHER':'##FISHER=<ID=*, Number=1, Type=Float, Description=\"Two-tailed Fisher exact test.\">\n',\
   'BOOST':'##BOOST=<ID=*, Number=1, Type=Float, Description=\"Scarcity-values, calculated from the allele frecuency values of the recessive bulk.\">\n',\
   'PVALUE':'##PVALUE=<ID=*, Number=1, Type=Float, Description=\"P-values of Fisher test.\">\n',\
   'log10PVALUE':'##log10PVALUE=<ID=*, Number=1, Type=Float, Description=\"Base 10 logarithm of p-value.\">\n',\
   'DELTA':'##DELTA=<ID=*, Number=1, Type=Float, Description=\"Difference between SNPidx1 minus SNPidx2 (1 to -1).\">\n',\
   'ED':'##ED=<ID=*, Number=1, Type=Float, Description=\"Euclidean distance.\">\n',\
   'G':'##G=<ID=*, Number=1, Type=Float, Description=\"G-statistic.\">\n',\
   'TYPE':'##TYPE=<ID=*, Number=0, Type=Flag, Description=\"Type of element in which the mutation occurs.\">\n',\
   'ID':'##ID=<ID=*, Number=0, Type=Flag, Description=\"ID of mutated element.\">\n',\
   'PARENT':'##PARENT=<ID=*, Number=0, Type=Flag, Description=\"The gene or mRNA containig the mutated element.\">\n',\
   'STRAND':'##STRAND=<ID=*, Number=0, Type=Flag, Description=\"Genome strand cointaining the element (+, -, .).\">\n',\
   'CODON_change':'##CODON_change=<ID=*, Number=0, Type=Flag, Description=\"Target codon before and after mutation (AGA > AAA).\">\n',\
   'AA_change':'##AA_change=<ID=*, Number=0, Type=Flag, Description=\"Aminoacid before and after mutation (S > F).\">\n',\
   'INFO':'##INFO=<ID=variant_type, Number=0, Type=Flag, Description=\"Substitution, Insertion or Deletion.\">\n'\
          '##INFO=<ID=effect, Number=0, Type=Flag, Description=\"Type of aminoacid change effect in CDS (Nonsynonymous/Synonymous:missense/nonsense/nonstop | Insertion/Deletion:in-frame/frameshift), in 5\'UTR (new_ATG:in_frame/out_of_frame:truncated_protein/elongated_protein), affects to non-coding gene (non_coding:exonic/intronic) or in deletions (genes/exons_deleted:id1,id2).\">\n'\
          '##INFO=<ID=left, Number=0, Type=Flag, Description=\"Gene or Exon immediately left (towards lower POS) and distance on bp(element:distance).\">\n'\
          '##INFO=<ID=right, Number=0, Type=Flag, Description=\"Gene or Exon immediately right (towards higher POS) and distance on bp(element:distance).\">\n'\
          '##INFO=<ID=5_splice_site, Number=0, Type=Flag, Description=\"The mutation is 3 or less bp from the upstream splicing site (exon or intron _boundary:ID_CDS).\">\n'\
          '##INFO=<ID=3_splice_site, Number=0, Type=Flag, Description=\"The mutation is 3 or less bp from the downstream splicing site (exon or intron _boundary:ID_CDS).\">\n'\
          '##INFO=<ID=5_prime_deletion, Number=0, Type=Flag, Description=\"Deletion of the 5-prime end of a gene (ID_gen:deletion_size).\">\n'\
          '##INFO=<ID=3_prime_deletion, Number=0, Type=Flag, Description=\"Deletion of the 3-prime end of a gene (ID_gen:deletion_size).\">\n'\
          '##INFO=<ID=5_splice_site_deletion, Number=0, Type=Flag, Description=\"Deletion of the 5 prime intron splice site (intron_boundary:nextExonName).\">\n'\
          '##INFO=<ID=3_splice_site_deletion, Number=0, Type=Flag, Description=\"Deletion of the 3 prime intron splice site (intron_boundary:prevExonName).\">\n'\
   
}
