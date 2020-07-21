## Purpose: writes out only variants passing all filters
## Usage: python print_pass_only.py <variants file with 'filter' column> <output filename>

import sys

outf = open(sys.argv[2],'w')

with open(sys.argv[1],'r') as f1:
	for line in f1:
		tmp = line.strip().split('\t')
		if tmp[0] == 'id':
			idx = {col:index for index, col in enumerate(tmp)}
			#print '\t'.join(tmp)
			outf.write('\t'.join(tmp) + '\n')
		else:
			filt = tmp[idx['filter']].split('|')
			keys = [i.split('_')[0] for i in filt]
			vals = [i.split('_')[1] for i in filt]

			filtd = dict(zip(keys, vals))

			if not 'FAIL' in tmp[idx['filter']]:
				#print '\t'.join(tmp)
				outf.write('\t'.join(tmp) + '\n')

outf.close()
