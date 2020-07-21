'''
Script to take as input sample table (downloaded from Terra), sample map (Picard), pedigree file, 
and output a 7-column tab-separated file containing
<sample_id>, 
<sample gvcf google bucket path>, <father gvcf path>, <mother gvcf path>
<sample gvcf index path>, <Father gvcf index path>, <mother gvcf index path>

'''

import sys
from optparse import OptionParser


####################################################################################################
## handle arguments
####################################################################################################
parser = OptionParser()
parser.add_option('-i', '--input', dest='sample_table',help='sample table (Terra)')
parser.add_option('-m', '--smap', dest='sample_map',help='sample map (picard)')
parser.add_option('-p', '--ped', dest='ped', help='pedigree file')
(options, args) = parser.parse_args()

## check all arguments present
if (options.sample_table == None or options.sample_map == None or options.ped == None):
	print('\n' + '## ERROR: missing arguments' + '\n')
	parser.print_help()
	print('\n')
	sys.exit()

sample_table = options.sample_table
sample_map = options.sample_map
ped = options.ped

## parse ped file
pedd = {} # { id : {'fa': father_id, 'mo': mother_id} }
with open(ped, 'r') as pedf:
  for line in pedf:
    tmp = line.strip().split('\t')
    fid = tmp[0]
    sid = tmp[1]
    faid = tmp[2]
    moid = tmp[3]
    sex = tmp[4] # 1=male, 2=female, other=unknown
    aff = tmp[5] # 1=unaffect, 2=affect, -9=missing, 0=missing
    pedd[sid] = {'fa': faid, 'mo': moid}

## parse sample_map
pathd = {} # {id : path}
with open(sample_map, 'r') as smapf:
  for line in smapf:
    tmp = line.strip().split('\t')
    id = tmp[0]
    path = tmp[1]
    pathd[id] = path


## parse sample table
index_d = {}
with open(sample_table, 'r') as stf:
	for line in stf:
		tmp = line.strip().split('\t')
		if tmp[0] == 'entity:sample_id':
			idx = {col:index for index, col in enumerate(tmp)}
		else:
			id = tmp[idx['entity:sample_id']]
			index = tmp[idx['output_vcf_index']]
			
			index_d[id] = index


## parse sample table
with open(sample_table, 'r') as stf2:
	for line in stf2:
		tmp = line.strip().split('\t')
		if tmp[0] == 'entity:sample_id':
			idx = {col:index for index, col in enumerate(tmp)}
			#print('\t'.join(['sample_id', 'sample_gvcf', 'father_gvcf', 'mother_gvcf', 'sample_gvcf_index', 'father_gvcf_index', 'mother_gvcf_index']))
		else:
			id = tmp[idx['entity:sample_id']]

			if not ('fa' in id or 'mo' in id):

				pb_path = pathd[id]
				pb_idx = index_d[id]
				
				if pedd[id]['fa'] == '0':
					fa_path = '.'
					fa_idx = '.'
				else:
					fa_path = pathd[pedd[id]['fa']]
					fa_idx = index_d[pedd[id]['fa']]
				if pedd[id]['mo'] == '0':
					mo_path = '.'
					mo_idx = '.'
				else:
					mo_path = pathd[pedd[id]['mo']]
					mo_idx = index_d[pedd[id]['mo']]


				if not ('fa' in id or 'mo' in id):
					print('\t'.join([id, pb_path, fa_path, mo_path, pb_idx, fa_idx, mo_idx]))
