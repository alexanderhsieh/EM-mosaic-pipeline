import sys
from optparse import OptionParser

########################################
## Handle arguments
########################################

parser = OptionParser()
parser.add_option('-i', '--input', dest='input_file',help='input tab-separated variants file')
parser.add_option('-c', '--cohort_AF_file', dest='cohort_AF_file', help='cohort allele frequency file')
parser.add_option('-o', '--output', dest='output_file',help='output filename')
(options, args) = parser.parse_args()

## check all arguments present
if (options.input_file == None or options.cohort_AF_file == None or options.output_file == None):
	print('\n' + '## ERROR: missing arguments' + '\n')
	parser.print_help()
	sys.exit()

input_file = options.input_file
cohort_AF_file = options.cohort_AF_file
output_file = options.output_file

########################################
## create dictionary of cohort AF
########################################
cafd = {} # { <chr:pos:ref:alt>: [carrier count, allele frequency] }
with open(cohort_AF_file, 'r') as cf:
	for line in cf:
		tmp = line.strip().split('\t')
		if tmp[0] == 'var_id':
			idx = {col:index for index, col in enumerate(tmp)}
		else:
			var_id = tmp[idx['var_id']]
			n_carriers = tmp[idx['n_carriers']]
			freq = tmp[idx['cohort_allele_frequency']]

			cafd[var_id] = [n_carriers, freq]

########################################
## append n_carriers and cohort_AF columns to corresponding variants
########################################
outf = open(output_file, 'w')

with open(input_file, 'r') as f:
	for line in f:
		tmp = line.strip().split('\t')
		if tmp[0] == 'id':
			idx = {col:index for index, col in enumerate(tmp)}
			head = '\t'.join(tmp) + '\t' + '\t'.join(['n_carriers', 'cohort_AF'])
			
			outf.write(head + '\n')

		else:
			chr = tmp[idx['chr']] # use vcf formatting for chromosome (chr1 instead of 1)
			pos = tmp[idx['pos']]
			ref = tmp[idx['ref']]
			alt = tmp[idx['alt']]

			tmp_var_id = ':'.join([chr, pos, ref, alt])

			out_n_carriers, out_cohort_AF = '.', '.'

			if tmp_var_id in cafd:
				tmpout = cafd[tmp_var_id]
				out_n_carriers = tmpout[0]
				out_cohort_AF = tmpout[1]


			out = '\t'.join(tmp) + '\t' + '\t'.join([out_n_carriers, out_cohort_AF])
			outf.write(out + '\n')

outf.close()
