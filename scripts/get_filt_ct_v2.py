## Purpose: parses filter column and counts the number of variants failing each respective filter
## Usage: python get_filt_ct.py <variants file with 'filter' column> <output filename>

import sys

outf = open(sys.argv[2],'w')

total = 0
allpass = 0
d = {'PV4':0, 'SB':0, 'FDR':0, 'RR':0, 'VC':0, 'MAF':0, 'COD':0, 'MUC-HLA':0, 'CAF':0, 'IGV':0, 'dbSNP':0, 'OUT':0}

with open(sys.argv[1],'r') as f1:
	for line in f1:
		tmp = line.strip().split('\t')
		if tmp[0] == 'id':
			idx = {col:index for index, col in enumerate(tmp)}
		else:
			total += 1

			filt = tmp[idx['filter']].split('|')
			keys = [i.split('_')[0] for i in filt]
			vals = [i.split('_')[1] for i in filt]

			filtd = dict(zip(keys, vals))

			for k in filtd:
				if filtd[k] == 'FAIL':
					d[k] += 1

			if not 'FAIL' in tmp[idx['filter']]:
				allpass += 1

#print 'TOTAL: %d'%(total)
outf.write('TOTAL: %d'%(total) + '\n')
for i in d:
	#print '#  %s : %d'%(i, d[i])
	outf.write('#  %s : %d'%(i, d[i]) + '\n')
#print 'PASSING ALL FILTERS: %d'%(allpass)
outf.write('PASSING ALL FILTERS: %d'%(allpass)+'\n')

outf.close()
