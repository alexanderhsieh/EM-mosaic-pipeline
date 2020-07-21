## Purpose: take raw variants file and format as ANNOVAR input file (.avinput)
## 			chr, start, end, ref, alt, otherinfo
'''
Usage: prep_avinput.py [options]

Options:
  -h, --help            show this help message and exit
  -i input_file, --input=input_file
                        input tab-separated variants file
  -c config_file, --config=config_file
                        configuration file
  -o output_file, --output=output_file
'''
import sys
from optparse import OptionParser
import os

## handle arguments
parser = OptionParser()
parser.add_option('-i', '--input', dest='input_file',help='input tab-separated variants file')
parser.add_option('-o', '--output', dest='output_file',help='output filename')
(options, args) = parser.parse_args()

## check all arguments present
if (options.input_file == None or  options.output_file == None):
	print('\n' + '## ERROR: missing arguments' + '\n')
	parser.print_help()
	sys.exit()


input_file = options.input_file
output_file = options.output_file

## dictionary of common header line names (to avoid needing yaml config file)
cfg_head = {
  'id': ['id', 'iid', 'individual', 'indv', 'sample', 'subject', 'individual_id', 'blinded id'],
  'chr': ['chr', 'chromosome', 'chrom', '#chrom'],
  'pos': ['pos', 'position', 'location', 'start'],
  'ref': ['ref', 'reference', 'ref_allele', 'ref.allele'],
  'alt': ['alt', 'alternate', 'alternative', 'alt_allele', 'alt.allele']
}


vcfout = open(output_file, 'w')

with open(input_file,'r') as f1:
	for line in f1:
		tmp = line.strip().split('\t')
		if not tmp[0].startswith('##'):
			## first col of header is either (1) VCF-style starting with #CHROM or (2) some id field
			if tmp[0].startswith('#') or tmp[0].lower() in cfg_head['id']:
				idx = {col:index for index, col in enumerate(tmp)} # original header dict
				idxmap = {} # re-mapped header dict
				
				## write out VCF header
				vcfout.write("##fileformat=VCFv4.1"+'\n')
				head = '\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'])
				#head = '\t'.join(['#CHROM','POS','REF','ALT','QUAL','FILTER','INFO','FORMAT'])
				vcfout.write(head+'\n')
			else:
				id = tmp[idx['id']]
				chr = tmp[idx['chr']]
				pos = tmp[idx['pos']]
				ref = tmp[idx['ref']]
				alt = tmp[idx['alt']]

				qual, filt, info, format = '.', '.','.','GT:AD:DP:GQ:PL'

				outstring = '\t'.join([chr, pos,id,  ref, alt, qual, filt, info, format])
				#outstring = '\t'.join([chr, pos, ref, alt, qual, filt, info, format])
				vcfout.write(outstring+'\n')

vcfout.close()


