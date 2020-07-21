## generates a table of <chr:pos:ref:alt> <number of samples carrying variant allele>
## estimated on the basis of genotypes

import sys
from optparse import OptionParser
import subprocess
import os
import gzip
import io

####################################################################################################
## handle arguments
####################################################################################################
parser = OptionParser()
parser.add_option('-i', '--input', dest='input_file',help='input annotated denovos file')
#parser.add_option('-m', '--min_vaf', dest='min_var',help='minimum VAF to be considered a carrier')
parser.add_option('-n', '--cohort_size', dest='cohort_size',help='number of samples in cohort')
parser.add_option('-o', '--output', dest='output_file',help='output tab-separated variants file')
(options, args) = parser.parse_args()

## check all arguments present
#if (options.input_gvcf == None or options.min_vaf == None or options.output_file):
if (options.input_file == None or options.cohort_size == None or options.output_file == None):
	print('\n' + '## ERROR: missing arguments' + '\n')
	parser.print_help()
	print('\n')
	sys.exit()


input_file = options.input_file
cohort_size = options.cohort_size
output_file = options.output_file


####################################################################################################
## parse gvcf
####################################################################################################
outf = open(output_file, 'w')
header = '\t'.join(['var_id', 'n_carriers', 'cohort_allele_frequency'])
outf.write(header + '\n')

carriers = {}
#with gzip.open(input_gvcf, 'r') as f:
with open(input_file, 'r') as f:
	for line in f:
		tmp = line.strip().split('\t')
		if tmp[0] == 'id':
			idx = {col:index for index, col in enumerate(tmp)}
		else:
			chr, pos, ref, alt = tmp[idx['chr']], tmp[idx['pos']], tmp[idx['ref']], tmp[idx['alt']]

			key = ':'.join([chr, pos, ref, alt])

			if not key in carriers:
				carriers[key] = 0
			else:
				carriers[key] += 1

with open(input_file, 'r') as f2:
	for line in f2:
		tmp = line.strip().split('\t')
		if tmp[0] == 'id':
			idx = {col:index for index, col in enumerate(tmp)}
		else:
			chr, pos, ref, alt = tmp[idx['chr']], tmp[idx['pos']], tmp[idx['ref']], tmp[idx['alt']]

			key = ':'.join([chr, pos, ref, alt])

			n_carriers = carriers[key]

			freq = float(n_carriers)/float(cohort_size)

			outf.write(key + '\t' + str(n_carriers) + '\t' + str(freq) + '\n')


outf.close()

