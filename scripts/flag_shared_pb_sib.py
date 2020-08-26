import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-i', '--input', dest='inputf', help='input file')
parser.add_option('-p', '--ped', dest='pedf', help='ped file')
parser.add_option('-o', '--out', dest='outf', help='output file')
(options, args) = parser.parse_args()

if (options.inputf == None or options.pedf == None or options.outf == None):
	print('\n'+'## ERROR: missing arguments')
	parser.print_help()
	print('\n')
	sys.exit()


## parse ped file to identify proband-sibling pairs
## map sample id to a key (father_id:mother_id)
## map key to list of children
famd = {} # { sid : faid:moid }
pedd = {} # { faid:moid : [children] }
with open(options.pedf, 'r') as pf:
	for line in pf:
		tmp = line.strip().split('\t')
		famid = tmp[0]
		sid = tmp[1]
		faid = tmp[2]
		moid = tmp[3]
		sex = tmp[4]
		aff = tmp[5]

		# ignore parent lines
		if not (faid == 0 and moid == 0):
			
			key = ':'.join([faid, moid])
			
			famd[sid] = key # sample_id to key mapping
			
			if not key in pedd:
				pedd[key] = []

			pedd[key].append(sid) # append sample_id to key


## first pass create dictionary of variants
tot = 0
vard = {} # { chr:pos:ref:alt : [carrier ids] }
with open(options.inputf, 'r') as f:
	for line in f:
		tmp = line.strip().split('\t')
		if tmp[0] == 'id':
			idx = {col:index for index, col in enumerate(tmp)}
		else:
			id = tmp[idx['id']]
			chr, pos, ref, alt = tmp[idx['chr']], tmp[idx['pos']], tmp[idx['ref']], tmp[idx['alt']]

			var = ':'.join([chr, pos, ref, alt])

			if not var in vard:
				vard[var] = []

			vard[var].append(id)

			tot += 1

## second pass

output_file = open(options.outf, 'w')

share_ct = 0
with open(options.inputf, 'r') as f:
	for line in f:
		tmp = line.strip().split('\t')
		if tmp[0] == 'id':
			idx = {col:index for index, col in enumerate(tmp)}
			head = '\t'.join(tmp) + '\t' + 'shared_pb_sib'
			output_file.write(head+'\n')
		else:
			id = tmp[idx['id']]
			chr, pos, ref, alt = tmp[idx['chr']], tmp[idx['pos']], tmp[idx['ref']], tmp[idx['alt']]

			var = ':'.join([chr, pos, ref, alt])

			## check if shared by proband sibling
			share_flag = 'False'
			isec = []

			## fetch list of children associated with sample's family
			tmp_famid = famd[id]
			tmp_children = pedd[tmp_famid]

			## get all variant carriers
			carriers = vard[var]

			## get intersection between variant carriers and proband-sibling list
			isec = list(set(carriers).intersection(set(tmp_children)))

			if len(isec) == 2:
				'''
				print('## variant')
				print(var)
				print('## children')
				print(tmp_children)
				print('## carriers')
				print(carriers)
				print('## intersect')
				print(isec)
				sys.exit()
				'''
				share_flag = 'True'
				share_ct += 1

			output_file.write('\t'.join(tmp) + '\t' + share_flag + '\n')


output_file.close()

print('%d/%d variants shared by proband and sibling pair'%(share_ct, tot))




