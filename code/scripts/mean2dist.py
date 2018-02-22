# import argparse
# parser = argparse.ArgumentParser(prog='mean2dist.py', description='''
#     Take a genotype mean file in BIMBAM format and output an imputed genotype distribution file where the posterior probability is either 0 or 1 and the closest genotype is set to one and others are set to zero. Also the alternative allele and reference allele are fliped.

#     Note that fliped genotype dist format is
#     [snp-id] [ref-allele] [alt-allele] [posterior-prob-of-ref-ref] [posterior-prob-of-ref-alt]
# ''')
# parser.add_argument('--input', help = 'input genotype mean file in GZ format')
# parser.add_argument('--output', help = 'file name of output genotype distribution file in GZ format')
# args = parser.parse_args()



input_f = snakemake.input[0]
output_f = snakemake.output[0]

import pandas as pd
import numpy as np

geno = pd.read_table(input_f, sep = ' ', header = None, compression = 'gzip')
geno_info = geno.iloc[:, [0, 2, 1]]
geno_mean = geno.iloc[:, 3:]
(ns, ni) = geno_mean.shape
geno_dist = pd.DataFrame()
distance = 0
for i in range(ni): # iter over every individual
    geno_mean_i = round(geno_mean.iloc[:, i])
    distance += sum(abs(geno_mean_i - geno_mean.iloc[:, i]))
    geno_dist_i = pd.DataFrame(np.zeros((ns, 2)))
    geno_dist_i.ix[geno_mean_i == 0, 0] = 1
    geno_dist_i.ix[geno_mean_i == 1, 1] = 1
    geno_dist = pd.concat([geno_dist, geno_dist_i], axis = 1)
geno_out = pd.concat([geno_info, geno_dist], axis = 1)
print('average distance between posterior mean genotype and integerized genotype is {avg_distance}'.format(avg_distance = distance / ns / ni))
geno_out.to_csv(output_f, sep = ' ', header = False, index = False, compression = 'gzip')
