## Purpose: parse VEP CSQ output and append designated columns to original denovos file
'''
Usage: parse_and_append_vep_cols.py [options] > output file

Options:
  -h, --help            show this help message and exit
  -i input_file, --input=input_file
                        input tab-separated variants file
  -v vep_file, --vep=vep_file
                        VEP-annotated VCF
  -c cols, --col=cols
                    	comma-separated list of columns to parse and append
  -o output_file, --output=output_file
                    	output tab-separated variants file with additional columns
'''
import sys
from optparse import OptionParser
import subprocess
import os


## handle arguments
parser = OptionParser()
parser.add_option('-i', '--input', dest='input_file',help='input tab-separated variants file')
parser.add_option('-v', '--vep', dest='vep_file',help='VEP-annotated VCF')
parser.add_option('-c', '--col', dest='cols', help='comma-separated list of columns to parse and append')
parser.add_option('-o', '--output', dest='output_file',help='output tab-separated variants file')
(options, args) = parser.parse_args()

## check all arguments present
if (options.input_file == None or options.vep_file == None or options.cols == None or options.output_file == None):
	print('\n' + '## ERROR: missing arguments' + '\n')
	parser.print_help()
	print('\n')
	sys.exit()


input_file = options.input_file
vep_file = options.vep_file
cols = options.cols
output_file = options.output_file

## dictionary of common header line names (to avoid needing yaml config file)
cfg_head = {
  'id': ['id', 'iid', 'individual', 'indv', 'sample', 'subject', 'individual_id', 'blinded id'],
  'chr': ['chr', 'chromosome', 'chrom', '#chrom'],
  'pos': ['pos', 'position', 'location', 'start'],
  'ref': ['ref', 'reference', 'ref_allele', 'ref.allele'],
  'alt': ['alt', 'alternate', 'alternative', 'alt_allele', 'alt.allele']
}

## parse VEP file and store as dict

vepd = {} # { id:chr:pos : csq } ; csq= { colname : colval }
with open(vep_file, 'r') as vf:
	for line in vf:
		## parse CSQ format
		if line.startswith('##INFO'):
			csqfmt = line.split('Format: ')[1].split('">')[0].split('|')
		## iterate over non-INFO lines
		elif not line.startswith('##'):
			tmp = line.strip().split('\t')
			if tmp[0] == '#CHROM':
				idx = {col:index for index, col in enumerate(tmp)}
			else:
				id = tmp[idx['ID']]
				chr, pos, ref, alt = tmp[idx['#CHROM']], tmp[idx['POS']], tmp[idx['REF']], tmp[idx['ALT']]

				csqval = tmp[idx['INFO']].split('CSQ=')[1].split('|')

				csq = dict(zip(csqfmt, csqval))

				#print '\n'.join(sorted(csq.keys()))
				#print '\t'.join([csq[c] for c in cols])

				key = ':'.join([id, chr, pos])

				vepd[key] = csq


## read input file and append designated columns

## list of output columns
outcols = cols.split(',')

## read header line to get first column value of header
#head = subprocess.Popen("head -n 1 %s"%(input_file), shell=True, stdout=subprocess.PIPE).stdout.read()
head = subprocess.run("head -n 1 %s"%(input_file), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout
head_token = head.split('\t')[0]

## added 03/26/2020; handles alternative input file format (VCF instead of txt)
if input_file.endswith('.vcf'):
	head_token = '#CHROM'


## added 04/16/2020: write to output_file instead of print to stdout()
outf = open(output_file, 'w')

with open(input_file, 'r') as f1:
	for line in f1:
		tmp = line.strip().split('\t')
		if line.startswith('##'): ## handle VCF header lines
			#print '\t'.join(tmp)
			outf.write('\t'.join(tmp) + '\n')
		else:
			if tmp[0] == head_token or tmp[0].lower() in cfg_head['id']:
				idx = {col:index for index, col in enumerate(tmp)} # original header dict

				outf.write('\t'.join(tmp) + '\t' + '\t'.join(outcols) + '\n')
			
			else:
				print(line)
				id = tmp[idx['id']]
				chr = tmp[idx['chr']]
				pos = tmp[idx['pos']]
				ref = tmp[idx['ref']]
				alt = tmp[idx['alt']]

				key = ':'.join([id, chr, pos])

				## fetch CSQ dictionary for this variant
				tmpd = vepd[key]

				## get user-specificed column values
				to_append = [tmpd[c] for c in outcols]

				## replace missing annotations with '.'
				to_append_fixed = ['.' if len(w) == 0 else w for w in to_append]  

				#print '\t'.join(tmp) + '\t' + '\t'.join(to_append_fixed)
				outf.write('\t'.join(tmp) + '\t' + '\t'.join(to_append_fixed) + '\n')


outf.close()
