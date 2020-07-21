## Purpose: parses BEDtools intersection results
## 			appends columns for mappability score, segdup score, LCR yes/no
##			note: since bedfiles are already parsed by score, default values: '.' indicating not hit by any
##			note: since some variants hit multiple regions with different scores, take max
## Usage: python parse_bedtools_isec.py <variants file> <bed.isec.out.txt>
## Output: prints to stdout()
import sys

## Read in bedtools intersect file and create dictionaries of scores
mapd = {} # { chr:pos : mappability score} 
segd = {} # { chr:pos : segdup score }
lcrd = {} # { chr:pos : LCR yes/no }

with open(sys.argv[2],'r') as bedf:
	for line in bedf:
		tmp = line.strip().split('\t')
		
		chr = tmp[0].strip('chr')
		pos = tmp[1]

		tmpkey = ':'.join([chr, pos])


		fname = tmp[3]


		# parse score columns where available
		if 'map' in fname:
			score = float(tmp[7])
			if not tmpkey in mapd:
				mapd[tmpkey] = score
			else:
				if score > mapd[tmpkey]:
					mapd[tmpkey] = score
		elif 'seg' in fname:
			score = float(tmp[7])
			if not tmpkey in segd:
				segd[tmpkey] = score
			else:
				if score > segd[tmpkey]:
					segd[tmpkey] = score
		elif ('LCR' in fname):
			lcrd[tmpkey] = 'yes'


## Update variants file and add new columns
with open(sys.argv[1],'r') as f1:
	for line in f1:
		tmp = line.strip().split('\t')
		if tmp[0] == 'id':
			idx = {col:index for index, col in enumerate(tmp)}
			print('\t'.join(tmp) + '\t' + '\t'.join(['map_score', 'segdup_score', 'LCR_flag']))
		else:
			chr, pos = tmp[idx['chr']], tmp[idx['pos']]

			key = ':'.join([chr, pos])

			out_map, out_sd, out_lcr = '.', '.', 'no'

			if key in mapd:
				out_map = mapd[key]
			if key in segd:
				out_sd = segd[key]
			if key in lcrd:
				out_lcr = lcrd[key]

			print('\t'.join(tmp) + '\t' + '\t'.join(map(str, [out_map, out_sd, out_lcr])))
